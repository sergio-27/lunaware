"""Run heat1d 0.3.2 and write the reference values used by lunarenv validation case TH-X1."""

import json, numpy as np, heat1d, planets

np.float = float
np.int = int  # aliases removed in NumPy 1.24; heat1d 0.3.2 still uses them
moon = planets.Moon
out = {
    "code": "heat1d",
    "version": heat1d.__version__ if hasattr(heat1d, "__version__") else "0.3.2",
    "source": "https://github.com/phayne/heat1d (PyPI heat1d 0.3.2)",
    "generated": "2026-11-12",
    "note": "Default grid (m=10, n=5; 21 layers; finer grids return NaN in this version). "
    "NYEARSEQ raised from the default 1 to 20 so that the deep temperature is equilibrated.",
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
    for years in (1, 20):
        m = heat1d.Model(
            planet=moon, lat=np.radians(lat), ndays=1, config=heat1d.Configurator(NYEARSEQ=years)
        )
        m.run()
        ts = m.T[:, 0]
        out["cases"].append(
            {
                "latitude_deg": lat,
                "equilibration_years": years,
                "sun_distance_au": float(m.r),
                "declination_rad": float(m.dec),
                "surface_mean_k": float(ts.mean()),
                "surface_max_k": float(ts.max()),
                "surface_min_k": float(ts.min()),
                "bottom_depth_m": float(m.profile.z[-1]),
                "bottom_mean_k": float(m.T[:, -1].mean()),
                "layers": int(m.N_z),
                "first_layer_m": float(m.profile.z[1]),
            }
        )
        print(out["cases"][-1])
json.dump(out, open("heat1d_reference.json", "w"), indent=1)
