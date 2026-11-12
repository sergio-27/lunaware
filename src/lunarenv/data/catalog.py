"""Catalog of external data sources and their pinned checksums.

The catalog (``catalog.toml``) says *what* can be downloaded and from where. The pin file
(``checksums.json``) says *exactly which bytes* are accepted. Keeping them apart lets the
pin file be machine-written while the catalog stays hand-curated.
"""

from __future__ import annotations

import json
import tomllib
from enum import StrEnum
from functools import cache
from importlib import resources
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, field_validator

try:  # advisory locking is POSIX-only
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None  # type: ignore[assignment]

_PACKAGE = "lunarenv.data"
_CATALOG_FILE = "catalog.toml"
_PINS_FILE = "checksums.json"


class DataKind(StrEnum):
    """File type; determines the signature check applied after download."""

    GEOTIFF = "geotiff"
    SPICE_LSK = "spice-lsk"
    SPICE_PCK_TEXT = "spice-pck-text"
    SPICE_PCK_BINARY = "spice-pck-binary"
    SPICE_FK = "spice-fk"
    SPICE_SPK = "spice-spk"


class CatalogEntry(BaseModel):
    """One downloadable file."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str = Field(pattern=r"^[a-z0-9_.]+$")
    url: str
    kind: DataKind
    groups: tuple[str, ...]
    product: str
    version: str
    citation: str
    description: str

    @field_validator("url")
    @classmethod
    def _https_only(cls, value: str) -> str:
        if not value.startswith("https://"):
            raise ValueError("data sources must use https")
        return value

    @property
    def relative_path(self) -> PurePosixPath:
        """Location inside the cache: ``<host>/<url path>`` (mirrors the source)."""
        parts = urlsplit(self.url)
        return PurePosixPath(parts.netloc) / parts.path.lstrip("/")


class Pin(BaseModel):
    """Accepted content of a catalog file."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    bytes: int = Field(gt=0)


class Catalog:
    """Immutable collection of catalog entries."""

    def __init__(self, entries: list[CatalogEntry]) -> None:
        self._entries: dict[str, CatalogEntry] = {}
        for entry in entries:
            if entry.id in self._entries:
                raise ValueError(f"duplicate catalog id: {entry.id}")
            self._entries[entry.id] = entry

    def __contains__(self, entry_id: str) -> bool:
        return entry_id in self._entries

    def __len__(self) -> int:
        return len(self._entries)

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(self._entries)

    def get(self, entry_id: str) -> CatalogEntry:
        try:
            return self._entries[entry_id]
        except KeyError:
            raise KeyError(f"unknown data id {entry_id!r}") from None

    def group(self, name: str) -> tuple[CatalogEntry, ...]:
        found = tuple(e for e in self._entries.values() if name in e.groups)
        if not found:
            raise KeyError(f"unknown data group {name!r}; known groups: {sorted(self.groups)}")
        return found

    @property
    def groups(self) -> frozenset[str]:
        return frozenset(g for e in self._entries.values() for g in e.groups)


def _expand_series(spec: dict[str, Any]) -> list[CatalogEntry]:
    first, last = int(spec.pop("first")), int(spec.pop("last"))
    if first > last:
        raise ValueError("series: first must not exceed last")
    entries = []
    for n in range(first, last + 1):
        fields = {k: v.format(n=n) if isinstance(v, str) else v for k, v in spec.items()}
        entries.append(CatalogEntry.model_validate(fields))
    return entries


def parse_catalog(text: str) -> Catalog:
    """Build a catalog from TOML text."""
    raw = tomllib.loads(text)
    unknown = set(raw) - {"file", "series"}
    if unknown:
        raise ValueError(f"unknown catalog sections: {sorted(unknown)}")
    entries = [CatalogEntry.model_validate(item) for item in raw.get("file", [])]
    for spec in raw.get("series", []):
        entries.extend(_expand_series(dict(spec)))
    return Catalog(entries)


@cache
def load_catalog() -> Catalog:
    """Return the catalog shipped with the package."""
    text = resources.files(_PACKAGE).joinpath(_CATALOG_FILE).read_text(encoding="utf-8")
    return parse_catalog(text)


def pins_path() -> Path:
    """Filesystem path of the pin file (writable only in a source checkout)."""
    return Path(str(resources.files(_PACKAGE).joinpath(_PINS_FILE)))


def load_pins(path: Path | None = None) -> dict[str, Pin]:
    """Read the pinned checksums. A missing file means nothing is pinned yet."""
    path = path or pins_path()
    if not path.exists():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {key: Pin.model_validate(value) for key, value in raw.items()}


def save_pin(entry_id: str, pin: Pin, path: Path | None = None) -> None:
    """Record a new pin. Refuses to change an existing one: pins are append-only.

    The read-modify-write is serialized across processes with an advisory lock (POSIX).
    """
    path = path or pins_path()
    with path.with_suffix(".lock").open("w") as lock:
        if fcntl is not None:
            fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            pins = load_pins(path)
            existing = pins.get(entry_id)
            if existing is not None and existing != pin:
                raise ValueError(
                    f"{entry_id} is already pinned to a different checksum; "
                    "edit checksums.json by hand if the upstream file was deliberately replaced"
                )
            pins[entry_id] = pin
            payload = {key: pins[key].model_dump() for key in sorted(pins)}
            tmp = path.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            tmp.replace(path)
        finally:
            if fcntl is not None:
                fcntl.flock(lock, fcntl.LOCK_UN)
