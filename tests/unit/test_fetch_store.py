"""Download verification and cache behaviour, against a local HTTP server."""

from __future__ import annotations

import hashlib
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from lunarenv.data import store as store_module
from lunarenv.data.catalog import Catalog, CatalogEntry, DataKind, Pin, load_pins
from lunarenv.data.fetch import DataIntegrityError, DownloadError, download
from lunarenv.data.store import DataNotAvailableError, DataStore, UnpinnedDataError

TIFF = b"II*\x00" + bytes(range(256)) * 8
HTML = b"<!DOCTYPE html><html><body>Not Found</body></html>"
BIG = b"II*\x00" + bytes(i * 7919 % 251 for i in range(20 << 20))
PAYLOADS = {"/good.tif": TIFF, "/missing.tif": HTML, "/empty.tif": b"", "/big.tif": BIG}


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/short.tif":
            self.send_response(200)
            self.send_header("Content-Length", str(len(TIFF) + 100))
            self.end_headers()
            self.wfile.write(TIFF)
            return
        body = PAYLOADS.get(self.path)
        if body is None:
            self.send_error(404)
            return
        byte_range = self.headers.get("Range")
        if byte_range and self.path == "/big.tif":
            first, last = (int(v) for v in byte_range.removeprefix("bytes=").split("-"))
            part = body[first : last + 1]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {first}-{last}/{len(body)}")
            self.send_header("Content-Length", str(len(part)))
            self.send_header("Last-Modified", "Fri, 05 Jun 2020 18:40:02 GMT")
            self.end_headers()
            self.wfile.write(part)
            return
        self.send_response(200)  # including for the HTML "not found" page, as PGDA does
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Last-Modified", "Fri, 05 Jun 2020 18:40:02 GMT")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture(scope="module")
def server():
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown()


def pin_of(data: bytes) -> Pin:
    return Pin(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))


def test_download_verifies_and_is_atomic(server, tmp_path):
    dest = tmp_path / "sub" / "good.tif"
    result = download(f"{server}/good.tif", dest, kind=DataKind.GEOTIFF, pin=pin_of(TIFF))
    assert dest.read_bytes() == TIFF
    assert result.sha256 == pin_of(TIFF).sha256 and result.bytes == len(TIFF)
    assert result.http_last_modified == "Fri, 05 Jun 2020 18:40:02 GMT"
    assert not list(dest.parent.glob("*.part"))


def test_html_served_with_status_200_is_rejected(server, tmp_path):
    dest = tmp_path / "missing.tif"
    with pytest.raises(DataIntegrityError, match="HTML"):
        download(f"{server}/missing.tif", dest, kind=DataKind.GEOTIFF)
    assert not dest.exists() and not list(tmp_path.glob("*.part"))


def test_pin_mismatch_is_rejected_and_nothing_is_kept(server, tmp_path):
    dest = tmp_path / "good.tif"
    with pytest.raises(DataIntegrityError, match="pinned"):
        download(f"{server}/good.tif", dest, kind=DataKind.GEOTIFF, pin=pin_of(b"other"))
    assert not dest.exists() and not list(tmp_path.glob("*.part"))


def test_truncated_and_absent_files_fail(server, tmp_path, monkeypatch):
    monkeypatch.setattr("time.sleep", lambda _: None)
    with pytest.raises(DownloadError):
        download(f"{server}/short.tif", tmp_path / "s.tif", kind=DataKind.GEOTIFF, retries=2)
    with pytest.raises(DownloadError):
        download(f"{server}/nope.tif", tmp_path / "n.tif", kind=DataKind.GEOTIFF, retries=2)
    with pytest.raises(DataIntegrityError):
        download(f"{server}/empty.tif", tmp_path / "e.tif", kind=DataKind.GEOTIFF)
    assert not list(tmp_path.iterdir())


def test_ranged_download_assembles_the_same_bytes(server, tmp_path):
    dest = tmp_path / "big.tif"
    result = download(
        f"{server}/big.tif", dest, kind=DataKind.GEOTIFF, pin=pin_of(BIG), connections=4
    )
    assert dest.read_bytes() == BIG and result.bytes == len(BIG)
    assert result.sha256 == pin_of(BIG).sha256


def test_ranged_download_falls_back_when_ranges_are_not_supported(server, tmp_path):
    dest = tmp_path / "good.tif"  # this path ignores Range and answers 200
    download(f"{server}/good.tif", dest, kind=DataKind.GEOTIFF, pin=pin_of(TIFF), connections=4)
    assert dest.read_bytes() == TIFF
    with pytest.raises(ValueError):
        download(f"{server}/good.tif", dest, kind=DataKind.GEOTIFF, connections=9)


# ------------------------------------------------------------------------- DataStore
@pytest.fixture
def local_store(server, tmp_path):
    entry = CatalogEntry.model_construct(  # bypass the https rule for the local server
        id="t.good",
        url=f"{server}/good.tif",
        kind=DataKind.GEOTIFF,
        groups=("g",),
        product="p",
        version="v",
        citation="c",
        description="d",
    )
    return DataStore(tmp_path / "cache", catalog=Catalog([entry]), pins_file=tmp_path / "pins.json")


def test_unpinned_files_are_refused_unless_pinning(local_store):
    with pytest.raises(UnpinnedDataError):
        local_store.fetch("t.good")
    with pytest.raises(DataNotAvailableError):
        local_store.path("t.good")
    record = local_store.fetch("t.good", pin=True)
    assert record.sha256 == pin_of(TIFF).sha256
    assert load_pins(local_store._pins_file)["t.good"] == pin_of(TIFF)
    assert local_store.path("t.good").read_bytes() == TIFF


def test_cached_file_is_returned_without_network(local_store, monkeypatch):
    first = local_store.fetch("t.good", pin=True)

    def no_network(*args, **kwargs):
        raise AssertionError("network used for a cached file")

    monkeypatch.setattr(store_module, "download", no_network)
    assert local_store.fetch("t.good") == first
    reopened = DataStore(
        local_store.root, catalog=local_store.catalog, pins_file=local_store._pins_file
    )
    assert reopened.fetch("t.good") == first  # the manifest persists


def test_present_file_without_record_is_adopted_if_it_matches_the_pin(local_store, monkeypatch):
    local_store.fetch("t.good", pin=True)
    path = local_store.path("t.good")
    local_store.manifest_path.unlink()
    fresh = DataStore(
        local_store.root, catalog=local_store.catalog, pins_file=local_store._pins_file
    )
    monkeypatch.setattr(store_module, "download", lambda *a, **k: pytest.fail("downloaded"))
    assert fresh.fetch("t.good").sha256 == pin_of(TIFF).sha256
    assert fresh.path("t.good") == path


def test_verify_detects_tampering_and_deletion(local_store):
    local_store.fetch("t.good", pin=True)
    assert local_store.verify() == {"t.good": "ok"}
    path = local_store.path("t.good")
    path.write_bytes(TIFF[:-1] + b"\x00")
    assert "differs" in local_store.verify()["t.good"]
    path.unlink()
    assert local_store.verify() == {"t.good": "missing"}


def test_tampered_cache_is_redownloaded_on_fetch(local_store):
    """A cached file whose record no longer matches the pin is replaced, not trusted."""
    local_store.fetch("t.good", pin=True)
    local_store.path("t.good").write_bytes(b"II*\x00corrupt")
    local_store.manifest_path.unlink()
    fresh = DataStore(
        local_store.root, catalog=local_store.catalog, pins_file=local_store._pins_file
    )
    fresh.fetch("t.good")
    assert fresh.path("t.good").read_bytes() == TIFF


def test_derived_products_keep_parent_checksum(local_store):
    parent = local_store.fetch("t.good", pin=True)
    out = local_store.new_derived_path("x/crop.bin")
    out.write_bytes(b"crop")
    record = local_store.register_derived(
        "t.good.crop",
        out,
        parent_id="t.good",
        operation="crop",
        parameters={"x_min": 0.0},
        created_utc="2026-01-01T00:00:00Z",
    )
    assert record.parent_sha256 == parent.sha256
    assert record.sha256 == hashlib.sha256(b"crop").hexdigest()
    assert local_store.derived_path("t.good.crop") == out
    local_store.discard_raw("t.good")
    assert not local_store.is_available("t.good")
    assert local_store.derived_path("t.good.crop").exists()  # provenance survives


def test_default_data_dir_order(tmp_path, monkeypatch):
    """Environment variable first, then the checkout's .data (as the tests), then the user cache."""
    monkeypatch.setenv(store_module.ENV_DATA_DIR, str(tmp_path / "env"))
    assert store_module.default_data_dir() == tmp_path / "env"
    monkeypatch.delenv(store_module.ENV_DATA_DIR)
    monkeypatch.setattr(store_module, "_checkout_data_dir", lambda: tmp_path / "checkout")
    assert store_module.default_data_dir() == tmp_path / "checkout"
    monkeypatch.setattr(store_module, "_checkout_data_dir", lambda: None)
    assert store_module.default_data_dir() == store_module.Path(
        store_module.user_cache_dir("lunarenv")
    )
