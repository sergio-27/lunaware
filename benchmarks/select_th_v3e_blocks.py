"""Select the TH-V3e Diviner blocks from model geometry (pre-registered 2026-11-17).

For each TH-V3e site, picks the 3 x 3 block of 240 m Diviner cells with the highest
summer mean visible solar disc at the cell centres (``diviner_pcp.select_block``) and
writes the selection to ``src/lunarenv/validation/reference/diviner_pcp_sites.json``
before any Diviner value is extracted.

Run:  LUNARENV_DATA_DIR=... uv run python benchmarks/select_th_v3e_blocks.py
"""

from __future__ import annotations

import json
import time
from datetime import UTC, datetime
from pathlib import Path

from lunarenv.data.store import DataStore
from lunarenv.validation import diviner_pcp as dp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src/lunarenv/validation/reference/diviner_pcp_sites.json"


def main() -> int:
    store = DataStore()
    selection = {}
    for label, near_id in dp.TH_V3E_SITES:
        tick = time.perf_counter()
        block = dp.select_block(store, near_id)
        selection[label] = block
        print(
            f"{label}: centre cell {block['centre_cell']}, latitude {block['latitude_deg']:.3f}, "
            f"visible {block['block_mean_visible_fraction']:.3f} ({time.perf_counter() - tick:.0f} s)",
            flush=True,
        )
    out = {
        "selected_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "rule": "docs/review/preregistration-th-v3e-2026-11-17.md, block selection",
        "selection": selection,
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
