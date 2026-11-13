"""Local data cache with a manifest of everything it contains.

The store is the only component that touches the network. Model code asks it for a path
by catalog id and receives a file whose checksum is known; results then embed the
corresponding manifest records as their data provenance.

Manifest updates are serialized across processes with an advisory file lock (POSIX), so
that parallel downloads do not lose each other's records.
"""

from __future__ import annotations

import os
import threading
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

from platformdirs import user_cache_dir
from pydantic import BaseModel, ConfigDict, Field

from lunarenv.data.catalog import (
    Catalog,
    CatalogEntry,
    DataKind,
    Pin,
    load_catalog,
    load_pins,
    save_pin,
)
from lunarenv.data.fetch import (
    MAX_CONNECTIONS,
    DataIntegrityError,
    check_signature,
    download,
    sha256_file,
    verify_pin,
)

try:  # advisory locking is POSIX-only; elsewhere updates are simply not serialized
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None  # type: ignore[assignment]

ENV_DATA_DIR = "LUNARENV_DATA_DIR"
_MANIFEST_FILE = "manifest.json"
_LOCK_FILE = ".manifest.lock"
_bookkeeping_lock = threading.RLock()  # pins and manifest, within one process


class DataNotAvailableError(RuntimeError):
    """A required file is not in the cache."""


class UnpinnedDataError(RuntimeError):
    """A file has no pinned checksum and pinning was not requested."""


class DataRecord(BaseModel):
    """Provenance of one cached source file."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str
    url: str
    relative_path: str
    kind: DataKind
    sha256: str
    bytes: int
    product: str
    version: str
    citation: str
    downloaded_utc: str
    http_last_modified: str | None = None


class DerivedRecord(BaseModel):
    """Provenance of a file derived from a source file (for example a crop)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str
    relative_path: str
    sha256: str
    bytes: int
    parent_id: str
    parent_sha256: str
    operation: str
    parameters: dict[str, float | int | str]
    created_utc: str


class Manifest(BaseModel):
    """Contents of a data cache."""

    model_config = ConfigDict(extra="forbid")

    schema_version: int = 1
    records: dict[str, DataRecord] = Field(default_factory=dict)
    derived: dict[str, DerivedRecord] = Field(default_factory=dict)


def _checkout_data_dir() -> Path | None:
    """``.data`` at the root of a source checkout, if the package runs from one and it exists."""
    root = Path(__file__).resolve().parents[3]  # src/lunarenv/data/store.py -> repository
    candidate = root / ".data"
    if (root / "pyproject.toml").is_file() and candidate.is_dir():
        return candidate
    return None


def default_data_dir() -> Path:
    """Cache location, in order: ``$LUNARENV_DATA_DIR``; ``.data`` at the root of the
    source checkout the package runs from, if that directory exists (the location the
    test suite uses); the per-user cache directory."""
    override = os.environ.get(ENV_DATA_DIR)
    if override:
        return Path(override).expanduser()
    return _checkout_data_dir() or Path(user_cache_dir("lunarenv"))


class DataStore:
    """Checksum-verified cache of the files declared in the catalog."""

    def __init__(
        self,
        root: Path | str | None = None,
        *,
        catalog: Catalog | None = None,
        pins_file: Path | None = None,
    ) -> None:
        self.root = (Path(root) if root is not None else default_data_dir()).resolve()
        self.catalog = catalog or load_catalog()
        self._pins_file = pins_file
        self._manifest = self._load_manifest()

    # ------------------------------------------------------------------ manifest
    @property
    def manifest_path(self) -> Path:
        return self.root / _MANIFEST_FILE

    @property
    def manifest(self) -> Manifest:
        return self._manifest

    def _load_manifest(self) -> Manifest:
        if self.manifest_path.exists():
            return Manifest.model_validate_json(self.manifest_path.read_text(encoding="utf-8"))
        return Manifest()

    @contextmanager
    def _manifest_lock(self) -> Iterator[None]:
        self.root.mkdir(parents=True, exist_ok=True)
        with (self.root / _LOCK_FILE).open("w") as handle:
            if fcntl is not None:
                fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                if fcntl is not None:
                    fcntl.flock(handle, fcntl.LOCK_UN)

    def _update_manifest(self, mutate: Callable[[Manifest], None]) -> None:
        """Apply ``mutate`` to the manifest on disk as one read-modify-write step."""
        with _bookkeeping_lock, self._manifest_lock():
            manifest = self._load_manifest()
            mutate(manifest)
            ordered = Manifest(
                records=dict(sorted(manifest.records.items())),
                derived=dict(sorted(manifest.derived.items())),
            )
            tmp = self.manifest_path.with_suffix(".json.tmp")
            tmp.write_text(ordered.model_dump_json(indent=2) + "\n", encoding="utf-8")
            tmp.replace(self.manifest_path)
            self._manifest = ordered

    # --------------------------------------------------------------------- paths
    def _raw_path(self, entry: CatalogEntry) -> Path:
        return self.root / "raw" / Path(*entry.relative_path.parts)

    def is_available(self, entry_id: str) -> bool:
        entry = self.catalog.get(entry_id)
        return entry_id in self._manifest.records and self._raw_path(entry).exists()

    def path(self, entry_id: str) -> Path:
        """Path of a cached file. Never downloads."""
        entry = self.catalog.get(entry_id)
        if not self.is_available(entry_id):
            raise DataNotAvailableError(
                f"{entry_id} is not in the data cache at {self.root}; "
                f"fetch it first (DataStore.fetch({entry_id!r}))"
            )
        return self._raw_path(entry)

    def record(self, entry_id: str) -> DataRecord:
        self.path(entry_id)
        return self._manifest.records[entry_id]

    def provenance(self, entry_ids: Iterable[str]) -> list[DataRecord]:
        """Manifest records for a set of ids, in a stable order."""
        return [self.record(entry_id) for entry_id in sorted(set(entry_ids))]

    # --------------------------------------------------------------------- fetch
    def fetch(
        self, entry_id: str, *, pin: bool = False, force: bool = False, connections: int = 1
    ) -> DataRecord:
        """Make a catalog file available locally and return its record.

        A file already in the cache whose recorded checksum matches the pin is returned
        without network access. With ``pin=True`` a file that has no pinned checksum is
        downloaded and its checksum recorded as the pin (maintainer operation, needs a
        writable source checkout). Without it, unpinned files are refused.
        """
        entry = self.catalog.get(entry_id)
        pins = load_pins(self._pins_file)
        pinned = pins.get(entry_id)
        if pinned is None and not pin:
            raise UnpinnedDataError(
                f"{entry_id} has no pinned checksum; a maintainer must fetch it with pin=True"
            )
        dest = self._raw_path(entry)
        if not force and self.is_available(entry_id):
            cached = self._manifest.records[entry_id]
            if pinned is None or (cached.sha256 == pinned.sha256 and cached.bytes == pinned.bytes):
                return cached
        if not force and pinned is not None and dest.exists():
            # File present without a record (interrupted run, copied cache). Adopt it if
            # its content is the pinned content; the download time is then unknown and
            # the file's modification time is recorded instead.
            sha256, size = sha256_file(dest), dest.stat().st_size
            if sha256 == pinned.sha256 and size == pinned.bytes:
                mtime = datetime.fromtimestamp(dest.stat().st_mtime, UTC)
                return self._store_record(
                    entry, sha256, size, mtime.strftime("%Y-%m-%dT%H:%M:%SZ"), None
                )
        result = download(entry.url, dest, kind=entry.kind, pin=pinned, connections=connections)
        if pinned is None:
            with _bookkeeping_lock:
                save_pin(entry_id, Pin(sha256=result.sha256, bytes=result.bytes), self._pins_file)
        return self._store_record(
            entry, result.sha256, result.bytes, result.downloaded_utc, result.http_last_modified
        )

    def _store_record(
        self,
        entry: CatalogEntry,
        sha256: str,
        size: int,
        downloaded_utc: str,
        http_last_modified: str | None,
    ) -> DataRecord:
        record = DataRecord(
            id=entry.id,
            url=entry.url,
            relative_path=str(entry.relative_path),
            kind=entry.kind,
            sha256=sha256,
            bytes=size,
            product=entry.product,
            version=entry.version,
            citation=entry.citation,
            downloaded_utc=downloaded_utc,
            http_last_modified=http_last_modified,
        )

        def mutate(manifest: Manifest) -> None:
            manifest.records[entry.id] = record

        self._update_manifest(mutate)
        return record

    def fetch_group(self, group: str, *, pin: bool = False) -> list[DataRecord]:
        """Fetch every catalog file that belongs to ``group``."""
        return [self.fetch(entry.id, pin=pin) for entry in self.catalog.group(group)]

    def fetch_many(
        self,
        entry_ids: Iterable[str],
        *,
        pin: bool = False,
        workers: int = 4,
        connections: int = 1,
    ) -> list[DataRecord]:
        """Fetch several files, ``workers`` at a time, each over ``connections`` streams.

        Transfers run in threads; each ``fetch`` call serializes its own pin and
        manifest update, so the bookkeeping stays consistent. The total number of
        simultaneous connections is capped: the sources are public science servers.
        """
        ids = list(dict.fromkeys(entry_ids))
        if workers < 1 or connections < 1 or workers * connections > MAX_CONNECTIONS:
            raise ValueError(f"workers x connections must be between 1 and {MAX_CONNECTIONS}")

        def one(entry_id: str) -> DataRecord:
            return self.fetch(entry_id, pin=pin, connections=connections)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(one, ids))

    def discard_raw(self, entry_id: str) -> None:
        """Delete a cached source file but keep its record (used after cropping)."""
        self._raw_path(self.catalog.get(entry_id)).unlink(missing_ok=True)

    # ------------------------------------------------------------------- derived
    def derived_path(self, derived_id: str) -> Path:
        """Path of a derived file. Raises if it is not registered."""
        try:
            record = self._manifest.derived[derived_id]
        except KeyError:
            raise DataNotAvailableError(
                f"no derived product {derived_id!r} in {self.root}"
            ) from None
        path = self.root / record.relative_path
        if not path.exists():
            raise DataNotAvailableError(f"derived product {derived_id!r} is registered but missing")
        return path

    def new_derived_path(self, relative_path: str) -> Path:
        """Reserve a location under ``derived/`` for a file about to be written."""
        path = self.root / "derived" / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def register_derived(
        self,
        derived_id: str,
        path: Path,
        *,
        parent_id: str,
        operation: str,
        parameters: dict[str, float | int | str],
        created_utc: str,
    ) -> DerivedRecord:
        """Record a derived file together with the checksum of its parent."""
        parent = self._manifest.records.get(parent_id)
        if parent is None:
            raise DataNotAvailableError(f"parent {parent_id!r} has no manifest record")
        record = DerivedRecord(
            id=derived_id,
            relative_path=str(path.resolve().relative_to(self.root)),
            sha256=sha256_file(path),
            bytes=path.stat().st_size,
            parent_id=parent_id,
            parent_sha256=parent.sha256,
            operation=operation,
            parameters=parameters,
            created_utc=created_utc,
        )

        def mutate(manifest: Manifest) -> None:
            manifest.derived[derived_id] = record

        self._update_manifest(mutate)
        return record

    # -------------------------------------------------------------------- verify
    def verify(self, entry_ids: Iterable[str] | None = None) -> dict[str, str]:
        """Re-hash cached files against manifest and pins.

        Returns a status per id: ``ok``, ``missing`` (record exists, file deleted) or a
        description of the mismatch. Nothing is raised, so a caller can report all
        problems at once.
        """
        pins = load_pins(self._pins_file)
        ids = sorted(self._manifest.records) if entry_ids is None else sorted(set(entry_ids))
        status: dict[str, str] = {}
        for entry_id in ids:
            record = self._manifest.records.get(entry_id)
            if record is None:
                status[entry_id] = "missing"
                continue
            path = self._raw_path(self.catalog.get(entry_id))
            if not path.exists():
                status[entry_id] = "missing"
                continue
            try:
                with path.open("rb") as handle:
                    check_signature(handle.read(16), record.kind)
                actual = Pin(sha256=sha256_file(path), bytes=path.stat().st_size)
                verify_pin(
                    record.sha256, record.bytes, actual, what=f"{entry_id} (manifest vs file)"
                )
                if entry_id in pins:
                    verify_pin(
                        actual.sha256,
                        actual.bytes,
                        pins[entry_id],
                        what=f"{entry_id} (file vs pin)",
                    )
            except DataIntegrityError as error:
                status[entry_id] = str(error)
            else:
                status[entry_id] = "ok"
        return status
