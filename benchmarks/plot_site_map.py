"""Figures for the illumination model documentation.

Reads .data/site_map_2026.npz (written by benchmarks/site_map.py) and writes
docs/figures/site_map_2026.png and docs/figures/site_horizon.png.

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/plot_site_map.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_series
from lunarenv.illumination.settings import IlluminationSettings, load_site_stack

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "docs" / "figures"


def main() -> int:
    store = DataStore()
    name = sys.argv[1] if len(sys.argv) > 1 else "site_map_2026"
    raw = np.load(store.root / f"{name}.npz")
    if "X" in raw:  # tile overview written by the coarse survey (100 m lattice)
        data = {
            "x": raw["X"][0],
            "y": raw["Y"][:, 0],
            "height": raw["h"],
            "mean": raw["mean"],
            "lcsp_d": raw["lcsp"] / 86400.0,
        }
        label = "Site 4 tile overview, one observer every 100 m"
    else:
        data = {k: raw[k] for k in ("x", "y", "height", "mean", "lcsp_d")}
        label = "Shackleton rim window, 5 m/pix"
    x, y = data["x"] / 1e3, data["y"] / 1e3
    extent = (x[0], x[-1], y[-1], y[0])
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    panels = [
        (data["height"], "Height above 1737.4 km sphere (m)", "terrain", None),
        (100 * data["mean"], "Mean visible solar disc fraction, 2026 (%)", "viridis", (0, 100)),
        (data["lcsp_d"], "Longest continuous shadow, 2026 (days)", "magma_r", None),
    ]
    for axis, (values, title, cmap, limits) in zip(axes, panels, strict=True):
        image = axis.imshow(
            values,
            extent=extent,
            origin="upper",
            cmap=cmap,
            vmin=None if limits is None else limits[0],
            vmax=None if limits is None else limits[1],
        )
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("polar stereographic x (km)")
        axis.set_ylabel("polar stereographic y (km)")
        fig.colorbar(image, ax=axis, shrink=0.85)
    best = np.unravel_index(int(np.argmax(data["mean"])), data["mean"].shape)
    for axis in axes:
        axis.plot(x[best[1]], y[best[0]], marker="+", color="red", markersize=12)
    fig.suptitle(f"{label}; observer 2 m above the surface", fontsize=11)
    fig.savefig(FIGURES / f"{name}.png", dpi=130)
    plt.close(fig)

    # horizon of the best pixel with the envelope of the Sun's path during 2026
    settings = IlluminationSettings()
    stack = load_site_stack(store, settings)
    bx, by = float(data["x"][best[1]]), float(data["y"][best[0]])
    horizons = compute_horizons(
        stack, bx, by, height_above_surface_m=2.0, n_azimuth=settings.n_azimuth
    )
    with spice_session(store) as ephemeris:
        utc = np.arange(
            np.datetime64("2026-01-01T00:00"),
            np.datetime64("2027-01-01T00:00"),
            np.timedelta64(60, "m"),
        )
        et = ephemeris.utc_to_et(utc)
        series = illumination_series(
            horizons, et, ephemeris.position_from_center("SUN", et), ephemeris.body_radius_m("SUN")
        )
    azimuth = np.degrees(horizons.azimuth_rad)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, constrained_layout=True)
    top.fill_between(
        azimuth, -3, np.degrees(horizons.elevation_rad[0]), color="0.35", label="terrain"
    )
    top.scatter(
        np.degrees(series.azimuth_rad),
        np.degrees(series.elevation_rad),
        s=1,
        c=series.visible_fraction,
        cmap="autumn",
        vmin=0,
        vmax=1,
        label="Sun centre, hourly (colour: visible fraction)",
    )
    top.set_ylim(-2.5, 3.0)
    top.set_ylabel("elevation (deg)")
    top.legend(loc="upper right", fontsize=8)
    top.set_title(
        f"Horizon at ({bx:.1f}, {by:.1f}) m, 2 m above the surface, and the Sun during 2026",
        fontsize=10,
    )
    bottom.semilogy(azimuth, horizons.distance_m[0] / 1e3, color="tab:blue", lw=0.8)
    bottom.set_ylabel("distance to skyline point (km)")
    bottom.set_xlabel("grid azimuth (deg, from map +y towards +x)")
    bottom.set_xlim(0, 360)
    fig.savefig(FIGURES / f"{name}_best_horizon.png", dpi=130)
    plt.close(fig)
    print("best pixel", bx, by, "mean", float(data["mean"][best]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
