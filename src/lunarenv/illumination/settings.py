"""Numerical settings of the illumination model and the terrain they act on.

Defaults were chosen from the convergence study in ``docs/convergence/illumination.md``;
the study, not this file, is the justification.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from lunarenv.data.dem import load_dem
from lunarenv.data.store import DataStore
from lunarenv.illumination.disc import CHORD, EDDINGTON_LIMB_DARKENING, PROFILE
from lunarenv.illumination.terrain import TerrainStack, build_stack

#: Numerical error of results at the default settings, from the convergence study
#: (docs/convergence/illumination.md). Estimates extrapolated from measured trends, not
#: bounds; they are reported next to the terrain-ensemble uncertainty.
NUMERICAL_ERROR_ESTIMATE = {
    "mean_illumination_percentage_points": 0.12,
    "mean_illumination_bias": "towards more illumination",
    "horizon_rms_arcsec": 65.0,
    "period_length": "one time step at each end",
}

#: Catalog ids of the terrain used for the MVP test site.
SITE04_NEAR = "pgda78.site04.surf"
REGIONAL_80M = "pgda90.ldem_80s_80mpp"


class IlluminationSettings(BaseModel):
    """Numerical parameters of horizon tracing and disc integration."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    n_azimuth: int = Field(1440, ge=8, description="Horizon directions per observer.")
    step_fraction: float = Field(
        0.5, gt=0.0, le=2.0, description="Ray sampling step in pixels of the level in use."
    )
    near_radius_m: float = Field(
        8000.0, gt=0.0, description="Ground distance out to which the site DEM is used."
    )
    mid_radius_m: float = Field(
        100_000.0, gt=0.0, description="Ground distance out to which the regional DEM is used."
    )
    far_coarsening: int = Field(
        3, ge=1, description="Block-averaging factor of the regional DEM beyond mid_radius_m."
    )
    n_strips: int = Field(32, ge=2, description="Strips used to integrate the disc.")
    limb_darkening: float = Field(
        EDDINGTON_LIMB_DARKENING, ge=0.0, le=1.0, description="Linear limb-darkening coefficient."
    )
    disc_mode: Literal["profile", "chord"] = Field(
        "profile", description="Horizon across the disc: per-strip profile, or one chord."
    )

    @property
    def disc_mode_code(self) -> int:
        return PROFILE if self.disc_mode == "profile" else CHORD


def load_site_stack(
    store: DataStore,
    settings: IlluminationSettings | None = None,
    *,
    near_id: str = SITE04_NEAR,
    regional_id: str = REGIONAL_80M,
) -> TerrainStack:
    """Terrain stack of a site from cached catalog products."""
    settings = settings or IlluminationSettings()
    return build_stack(
        load_dem(store, near_id),
        load_dem(store, regional_id),
        near_radius_m=settings.near_radius_m,
        mid_radius_m=settings.mid_radius_m,
        step_fraction=settings.step_fraction,
        far_coarsening=settings.far_coarsening,
    )
