"""Run heat1d and write the reference values used by lunarenv validation cases TH-X1 / TH-X1b.

heat1d 0.3.2 (PyPI, obsolete) -> heat1d_reference.json (TH-X1, kept as history).
heat1d 0.4.x (GitHub, current) -> heat1d_reference_current.json (TH-X1b), e.g.
heat1d 0.4.x with a second argument "converged" -> heat1d_reference_converged.json (TH-X1c:
grids m=n=20 and m=n=40, 20-year equilibration), e.g.

    uv run --isolated --with "heat1d @ git+https://github.com/phayne/heat1d.git@5b73a6c#subdirectory=python" \
        python benchmarks/make_heat1d_reference.py
"""

import json
import sys

import numpy as np

import heat1d

CURRENT = not getattr(heat1d, "__version__", "0.3.2").startswith("0.3")
CONVERGED = CURRENT and len(sys.argv) > 2 and sys.argv[2] == "converged"
if CURRENT:
    from heat1d import planets
else:
    import planets

    np.float = float
    np.int = int  # aliases removed in NumPy 1.24; heat1d 0.3.2 still uses them
moon = planets.Moon
out = {
    "code": "heat1d",
    "version": heat1d.__version__ if hasattr(heat1d, "__version__") else "0.3.2",
    "source": (
        "https://github.com/phayne/heat1d (GitHub, commit given on the command line)"
        if CURRENT
        else "https://github.com/phayne/heat1d (PyPI heat1d 0.3.2)"
    ),
    "commit": sys.argv[1] if len(sys.argv) > 1 else None,
    "generated": "2026-11-16" if CONVERGED else "2026-11-15" if CURRENT else "2026-11-12",
    "note": (
        "Grids m=n=20 and m=n=40 (TH-X1c pre-registration, 2026-11-16), default solver, "
        "NYEARSEQ 20; the m=n=40 grid is the reference. Means time-weighted over local "
        "time; the plain sample mean is also stored."
        if CONVERGED
        else "Default grid and solver of this version. equilibration_years 20 = NYEARSEQ 20 (as "
        "TH-X1, the primary reference); 0 = the version's automatic equilibration (DTBOT "
        "convergence), kept for information. Choice fixed before the run (TH-X1b "
        "pre-registration). Means are time-weighted (trapezoid over local time) because "
        "this version's output steps are not uniform; the plain sample mean is also stored."
        if CURRENT
        else "Default grid (m=10, n=5; 21 layers; finer grids return NaN in this version). "
        "NYEARSEQ raised from the default 1 to 20 so that the deep temperature is equilibrated."
    ),
    "planet": {
        "albedo": moon.albedo,
        "albedo_coefficients": list(moon.albedoCoef),
        "emissivity": moon.emissivity,
        "ks": moon.ks,
        "kd": moon.kd,
        "rhos": moon.rhos,
        "rhod": moon.rhod,
        "H": moon.H,
        "Qb": moon.Qb,
        "day_s": moon.day,
        "S": moon.S,
    },
    "cases": [],
}
for lat in (0.0, 20.0, 26.0):
    if CONVERGED:
        variants = [{"NYEARSEQ": 20, "m": k, "n": k} for k in (20, 40)]
    else:
        variants = [{"NYEARSEQ": y} for y in ((0, 20) if CURRENT else (1, 20))]
    for config in variants:
        years = config["NYEARSEQ"]
        m = heat1d.Model(
            planet=moon, lat=np.radians(lat), ndays=1, config=heat1d.Configurator(**config)
        )
        m.run()
        ts = m.T[:, 0]
        if CURRENT:
            # 0.4.x writes output at non-uniform local-time steps (denser by day), so a plain
            # sample mean is biased warm: weight by time instead.
            lt = np.asarray(m.lt, dtype=float)
            span = lt[-1] - lt[0]
            surface_mean = float(np.trapezoid(ts, lt) / span)
            bottom_mean = float(np.trapezoid(m.T[:, -1], lt) / span)
        else:
            surface_mean = float(ts.mean())
            bottom_mean = float(m.T[:, -1].mean())
        out["cases"].append(
            {
                "latitude_deg": lat,
                "equilibration_years": years,
                "grid_m_n": [config["m"], config["n"]] if CONVERGED else None,
                "sun_distance_au": float(m.r),
                "declination_rad": float(m.dec),
                "surface_mean_k": surface_mean,
                "surface_mean_k_plain_sample_mean": float(ts.mean()),
                "surface_max_k": float(ts.max()),
                "surface_min_k": float(ts.min()),
                "bottom_depth_m": float(m.profile.z[-1]),
                "bottom_mean_k": bottom_mean,
                "layers": int(m.N_z),
                "first_layer_m": float(m.profile.z[1]),
            }
        )
        print(out["cases"][-1])
name = (
    "heat1d_reference_converged.json"
    if CONVERGED
    else "heat1d_reference_current.json"
    if CURRENT
    else "heat1d_reference.json"
)
json.dump(out, open(name, "w"), indent=1)
print("wrote", name)
