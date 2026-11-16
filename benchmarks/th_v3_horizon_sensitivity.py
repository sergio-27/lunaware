"""Sensitivity of TH-V3 / TH-V3b to the resolution of the ground-level horizon.

Informational, not a validation case. Runs the TH-V3 forward model (9 Diviner cells of
240 m, 20 m pixel spacing) with the horizons computed on the 5 m site DEM (as TH-V3) and
on that DEM block-averaged to 10 and 20 m; the surface normal stays on the 5 m DEM. For
each albedo law it reports the block bias over the sunlit summer bins, both with the
sunlit bins of each run and with the bins of the 5 m run, and the per-cell bias.

    LUNARENV_DATA_DIR=$PWD/.data uv run python benchmarks/th_v3_horizon_sensitivity.py

About 20 minutes per resolution on 8 worker processes. Output:
docs/convergence/th_v3_horizon_sensitivity.json
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from lunarenv.data.store import DataStore
from lunarenv.validation import diviner_pcp as dp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/convergence/th_v3_horizon_sensitivity.json"
FACTORS = (1, 2, 4)  # 5, 10, 20 m horizons


def _tbol(radiance: np.ndarray, count: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray((radiance / count) ** 0.25)


def main() -> int:
    store = DataStore(os.environ.get("LUNARENV_DATA_DIR", str(ROOT / ".data")))
    reference = dp.load_reference()
    keys = [c["key"] for c in reference["cells"]]
    centres = [(float(c["x_m"]), float(c["y_m"])) for c in reference["cells"]]
    diviner = np.full((len(keys), dp.N_BINS), np.nan)
    for b, entry in enumerate(reference["bins"]):
        for c, k in enumerate(keys):
            if k in entry["tbol_k"]:
                diviner[c, b] = entry["tbol_k"][k]
    complete = np.all(np.isfinite(diviner), axis=0)
    diviner_block = np.where(
        complete, np.mean(np.where(complete, diviner, 0.0) ** 4, axis=0) ** 0.25, np.nan
    )
    runs: dict[str, dict] = {}  # type: ignore[type-arg]
    baseline_used = None
    for factor in FACTORS:
        cells = dp.forward(store, centres, 20.0, horizon_coarsening=factor)
        visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)
        used = (visible >= 0.5) & complete
        if baseline_used is None:
            baseline_used = used
        entry: dict[str, object] = {
            "horizon_dem_m": 5 * factor,
            "sunlit_bins": int(used.sum()),
            "mean_visible_fraction_summer": float(np.mean(visible[complete])),
        }
        for law in dp.LAWS:
            model_cell = np.array([_tbol(c.sums[law], c.sums["count"]) for c in cells])
            model_block = np.mean(model_cell**4, axis=0) ** 0.25
            diff = model_block - diviner_block
            entry[law] = {
                "block_bias_k": float(np.mean(diff[used])),
                "block_bias_k_on_5m_bins": float(np.mean(diff[baseline_used])),
                "cell_bias_k": {
                    k: float(np.nanmean((model_cell[c] - diviner[c])[baseline_used]))
                    for c, k in enumerate(keys)
                },
            }
        runs[f"{5 * factor}m"] = entry
        print(json.dumps(entry)[:400], flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps({"note": __doc__.strip().splitlines()[0], "runs": runs}, indent=1) + "\n"
    )
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
