"""Download a file and verify what arrived.

Verification is layered, because each layer catches a different failure:

1. *Signature*: the first bytes must match the declared file type. Some servers answer a
   missing file with ``200 OK`` and an HTML page (observed on pgda.gsfc.nasa.gov), so the
   HTTP status alone is not evidence that the data arrived.
2. *Length*: bytes written must equal ``Content-Length`` when the server sends it.
3. *Checksum*: SHA-256 must equal the pinned value when one exists.

The file only appears at its final path after all checks pass (atomic rename).

Large files can be fetched over several connections with HTTP range requests when the
server supports them (some servers limit the rate per connection). The checks are the
same: the checksum is computed over the assembled file.
"""

from __future__ import annotations

import hashlib
import os
import ssl
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import truststore

from lunarenv import __version__
from lunarenv.data.catalog import DataKind, Pin

_CHUNK_BYTES = 1 << 20
_MIN_SEGMENT_BYTES = 8 << 20
#: Simultaneous connections per process. pgda.gsfc.nasa.gov limits each connection to
#: roughly 100 kB/s and answered 503 when eight files were requested at once.
MAX_CONNECTIONS = 8
_BUSY_STATUS = (429, 503)  # the server asks clients to slow down
_BUSY_WAIT_S = 30.0
_USER_AGENT = f"lunarenv/{__version__}"

_SIGNATURES: dict[DataKind, tuple[bytes, ...]] = {
    # classic TIFF and BigTIFF, little- and big-endian
    DataKind.GEOTIFF: (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"),
    DataKind.SPICE_LSK: (b"KPL/LSK",),
    DataKind.SPICE_PCK_TEXT: (b"KPL/PCK",),
    DataKind.SPICE_FK: (b"KPL/FK",),
    DataKind.SPICE_PCK_BINARY: (b"DAF/PCK",),
    DataKind.SPICE_SPK: (b"DAF/SPK",),
}


class DownloadError(RuntimeError):
    """The file could not be retrieved."""


class DataIntegrityError(RuntimeError):
    """The file was retrieved but is not what was expected."""


@dataclass(frozen=True, slots=True)
class DownloadResult:
    sha256: str
    bytes: int
    downloaded_utc: str
    http_last_modified: str | None


def _ssl_context() -> ssl.SSLContext:
    """TLS context that verifies against the operating system trust store.

    Standalone Python builds (such as those installed by uv) ship without root
    certificates, and organizations often add their own; the OS store covers both.
    """
    return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)


def sha256_file(path: Path) -> str:
    """SHA-256 of a file, streamed."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def check_signature(head: bytes, kind: DataKind) -> None:
    """Raise if ``head`` (the first bytes of a file) does not match ``kind``."""
    if not any(head.startswith(sig) for sig in _SIGNATURES[kind]):
        looks_html = head.lstrip()[:1] == b"<"
        hint = " (the server returned an HTML page: the file probably does not exist)"
        raise DataIntegrityError(
            f"content does not look like {kind.value}: first bytes {head[:8]!r}"
            + (hint if looks_html else "")
        )


def download(
    url: str,
    dest: Path,
    *,
    kind: DataKind,
    pin: Pin | None = None,
    retries: int = 3,
    timeout_s: float = 60.0,
    connections: int = 1,
) -> DownloadResult:
    """Download ``url`` to ``dest`` and verify it.

    Network failures are retried with exponential backoff. Integrity failures are not:
    a wrong file will be wrong again. With ``connections > 1`` the file is fetched in
    byte ranges in parallel if the server allows it, and in one stream otherwise.
    """
    if not 1 <= connections <= MAX_CONNECTIONS:
        raise ValueError(f"connections must be between 1 and {MAX_CONNECTIONS}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    last_error: Exception | None = None
    for attempt in range(retries):
        if attempt:
            time.sleep(2.0**attempt)
        try:
            result = None
            if connections > 1:
                result = _download_ranged(
                    url, partial, kind=kind, timeout_s=timeout_s, connections=connections
                )
            if result is None:
                result = _download_once(url, partial, kind=kind, timeout_s=timeout_s)
        except urllib.error.HTTPError as error:
            partial.unlink(missing_ok=True)
            if error.code in _BUSY_STATUS:
                time.sleep(_BUSY_WAIT_S * (attempt + 1))
            elif 400 <= error.code < 500:  # the request itself is wrong; retrying cannot help
                raise DownloadError(f"{url}: HTTP {error.code} {error.reason}") from error
            last_error = error
            continue
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            last_error = error
            partial.unlink(missing_ok=True)
            continue
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
        try:
            if pin is not None:
                verify_pin(result.sha256, result.bytes, pin, what=url)
        except DataIntegrityError:
            partial.unlink(missing_ok=True)
            raise
        partial.replace(dest)
        return result
    raise DownloadError(f"could not download {url} after {retries} attempts: {last_error}")


def verify_pin(sha256: str, size: int, pin: Pin, *, what: str) -> None:
    """Raise if a checksum/size pair does not match its pin."""
    if size != pin.bytes or sha256 != pin.sha256:
        raise DataIntegrityError(
            f"{what}: content differs from the pinned version "
            f"(got {size} bytes, sha256 {sha256}; expected {pin.bytes} bytes, sha256 {pin.sha256})"
        )


def _download_once(url: str, partial: Path, *, kind: DataKind, timeout_s: float) -> DownloadResult:
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    digest = hashlib.sha256()
    written = 0
    with urllib.request.urlopen(request, timeout=timeout_s, context=_ssl_context()) as response:
        declared = response.headers.get("Content-Length")
        last_modified = response.headers.get("Last-Modified")
        with partial.open("wb") as handle:
            first = True
            while chunk := response.read(_CHUNK_BYTES):
                if first:
                    check_signature(chunk[:16], kind)
                    first = False
                handle.write(chunk)
                digest.update(chunk)
                written += len(chunk)
    if written == 0:
        raise DataIntegrityError(f"{url}: empty response")
    if declared is not None and int(declared) != written:
        # a short read is a transport problem, so let the caller retry
        raise ConnectionError(f"{url}: received {written} of {declared} bytes")
    return DownloadResult(
        sha256=digest.hexdigest(),
        bytes=written,
        downloaded_utc=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        http_last_modified=last_modified,
    )


def _open(url: str, timeout_s: float, byte_range: tuple[int, int] | None = None):  # type: ignore[no-untyped-def]
    headers = {"User-Agent": _USER_AGENT}
    if byte_range is not None:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    request = urllib.request.Request(url, headers=headers)
    return urllib.request.urlopen(request, timeout=timeout_s, context=_ssl_context())


def _download_ranged(
    url: str, partial: Path, *, kind: DataKind, timeout_s: float, connections: int
) -> DownloadResult | None:
    """Fetch a file in parallel byte ranges. Returns ``None`` if ranges are not supported."""
    with _open(url, timeout_s, (0, 15)) as probe:
        content_range = probe.headers.get("Content-Range", "")
        if probe.status != 206 or "/" not in content_range:
            return None
        total = int(content_range.rsplit("/", 1)[1])
        last_modified = probe.headers.get("Last-Modified")
        check_signature(probe.read(16), kind)
    count = max(1, min(connections, total // _MIN_SEGMENT_BYTES))
    bounds = [(i * total // count, (i + 1) * total // count - 1) for i in range(count)]

    descriptor = os.open(partial, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        os.ftruncate(descriptor, total)

        def fetch_segment(segment: tuple[int, int]) -> None:
            first, last = segment
            error: Exception | None = None
            for attempt in range(3):
                if attempt:
                    time.sleep(2.0**attempt)
                try:
                    offset = first
                    with _open(url, timeout_s, (first, last)) as response:
                        if response.status != 206:
                            raise ConnectionError("server stopped honouring range requests")
                        while chunk := response.read(_CHUNK_BYTES):
                            os.pwrite(descriptor, chunk, offset)
                            offset += len(chunk)
                    if offset != last + 1:
                        raise ConnectionError(
                            f"segment ended at byte {offset}, expected {last + 1}"
                        )
                    return
                except (urllib.error.URLError, TimeoutError, ConnectionError) as caught:
                    error = caught
            raise ConnectionError(f"{url}: segment {first}-{last} failed: {error}")

        with ThreadPoolExecutor(max_workers=count) as pool:
            list(pool.map(fetch_segment, bounds))
    finally:
        os.close(descriptor)
    return DownloadResult(
        sha256=sha256_file(partial),
        bytes=total,
        downloaded_utc=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        http_last_modified=last_modified,
    )
