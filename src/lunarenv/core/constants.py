"""Physical and geodetic constants, each with its unit and source.

Constants are plain floats in SI units so that they can be passed to compiled kernels.
Their provenance lives in :data:`CONSTANTS`, which is exported with every model
description so that a reader of a result can see exactly which values were used.

Reference keys resolve in ``docs/references.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True, slots=True)
class Constant:
    """A named constant with unit and provenance."""

    name: str
    value: float
    unit: str
    source: str
    note: str = ""


#: Radius of the reference sphere of the LOLA polar elevation models. Heights in the
#: DEMs are radial distances from the Moon's centre of mass minus this value.
MOON_REFERENCE_RADIUS_M: Final = 1_737_400.0

#: Stefan-Boltzmann constant (exact in the 2019 SI).
STEFAN_BOLTZMANN_W_M2_K4: Final = 5.670374419e-8

#: Astronomical unit (exact, IAU 2012 Resolution B2).
ASTRONOMICAL_UNIT_M: Final = 149_597_870_700.0

#: Total solar irradiance at 1 au.
SOLAR_CONSTANT_W_M2: Final = 1361.0

#: Speed of light in vacuum (exact).
SPEED_OF_LIGHT_M_S: Final = 299_792_458.0


CONSTANTS: Final[dict[str, Constant]] = {
    c.name: c
    for c in (
        Constant(
            "moon_reference_radius",
            MOON_REFERENCE_RADIUS_M,
            "m",
            "PGDA-78; verified against the CRS of every DEM at read time",
            "Sphere; planetocentric and planetodetic latitude coincide.",
        ),
        Constant(
            "stefan_boltzmann",
            STEFAN_BOLTZMANN_W_M2_K4,
            "W m-2 K-4",
            "CODATA 2018 (exact in the 2019 SI)",
        ),
        Constant("astronomical_unit", ASTRONOMICAL_UNIT_M, "m", "IAU 2012 Resolution B2 (exact)"),
        Constant(
            "solar_constant",
            SOLAR_CONSTANT_W_M2,
            "W m-2",
            "Hayne2017 Table A1, citing KoppLean2011",
        ),
        Constant("speed_of_light", SPEED_OF_LIGHT_M_S, "m s-1", "SI definition (exact)"),
    )
}
