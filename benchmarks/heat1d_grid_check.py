"""Post-hoc check behind the TH-X1b investigation: sensitivity of heat1d 0.4.x to its solver and grid.

Equator, 20-year equilibration, lunar defaults; prints the time-weighted surface mean, the
maximum and the minimum for heat1d's default settings, other solvers and finer grids.

    uv run --isolated --python 3.12 --with \
        "heat1d @ git+https://github.com/phayne/heat1d.git@5b73a6c#subdirectory=python" \
        python benchmarks/heat1d_grid_check.py
"""

import numpy as np

import heat1d
from heat1d import planets


def run(**settings: object) -> None:
    m = heat1d.Model(
        planet=planets.Moon, lat=0.0, ndays=1, config=heat1d.Configurator(NYEARSEQ=20, **settings)
    )
    m.run()
    ts = m.T[:, 0]
    lt = np.asarray(m.lt, dtype=float)
    mean = np.trapezoid(ts, lt) / (lt[-1] - lt[0])
    print(
        settings or "default",
        f"layers {m.N_z}",
        f"mean {mean:.2f} max {ts.max():.2f} min {ts.min():.2f}",
    )


run()
run(adaptive_tol=None)
run(solver="crank-nicolson")
run(solver="implicit", adaptive_tol=None)
for k in (10, 20, 40):
    run(m=k, n=k)
