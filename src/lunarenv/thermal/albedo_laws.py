"""Named albedo laws and the polar albedo envelope (scientific review, change 2).

The albedo at grazing incidence is the largest model-form uncertainty of the polar
thermal result (review finding F03). Instead of one law, results are given for an
envelope of published highland laws, read from the primary sources on 2026-11-16:

* cold bound, ``burger2024``: Bürger et al. (2024), JGR Planets 129, e2023JE008152,
  Eq. 15 with the section 4.2.2 best fit for the highlands high-conductivity model,
  A0 = 0.09, a = 0.10, b = 0.02 (0.860 at 88.2 degrees). Fitted to Diviner up to 80
  degrees latitude; steepest rise at grazing incidence;
* middle, ``hayne2017``: Hayne et al. (2017) Table A1, the engine default (0.784);
* hot bound, ``foote2020``: Foote et al. (2020), Icarus 336, 113456, Eq. 9, Apollo 16
  highland soil, laboratory hemispheric albedo (0.238). Measured to 75 degrees
  incidence; extrapolated beyond.

The review attributes the same two bounds to the DSNE Rev J highland envelope, which
could not be obtained. The constant albedo 0.12 (``constant``) lies outside the envelope
on the hot side; it is kept as the separate alternative TH-V3b tested.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from lunarenv.thermal.properties import RegolithProperties


@dataclass(frozen=True, slots=True)
class AlbedoLaw:
    """An albedo law as a set of RegolithProperties fields, with its source."""

    name: str
    source: str
    note: str
    fields: dict[str, Any] = field(default_factory=dict)

    def apply(self, properties: RegolithProperties) -> RegolithProperties:
        """The same properties with this albedo law (all albedo fields replaced)."""
        update: dict[str, Any] = {
            "albedo_normal": properties.albedo_normal,
            "albedo_a": properties.albedo_a,
            "albedo_b": properties.albedo_b,
            "albedo_quadratic_deg": None,
        }
        update.update(self.fields)
        return properties.model_copy(update=update)


BURGER2024 = AlbedoLaw(
    "burger2024",
    "Bürger et al. (2024) JGR Planets 129, e2023JE008152, Eq. 15, section 4.2.2",
    "fit to Diviner temperatures up to 80 deg latitude, highlands; cold bound",
    {"albedo_normal": 0.09, "albedo_a": 0.10, "albedo_b": 0.02},
)
HAYNE2017 = AlbedoLaw(
    "hayne2017",
    "Hayne et al. (2017) JGR Planets 122, Table A1",
    "fit to equatorial Diviner temperatures; engine default",
    {"albedo_normal": 0.12, "albedo_a": 0.06, "albedo_b": 0.25},
)
FOOTE2020 = AlbedoLaw(
    "foote2020",
    "Foote et al. (2020) Icarus 336, 113456, Eq. 9 (Apollo 16 highland soil)",
    "laboratory hemispheric albedo, measured to 75 deg incidence and extrapolated beyond; "
    "hot bound",
    {"albedo_quadratic_deg": (0.14382, -1.0248e-4, 1.3299e-5)},
)
CONSTANT = AlbedoLaw(
    "constant",
    "Hayne et al. (2017) A0 without incidence dependence",
    "outside the envelope on the hot side; the alternative of TH-V3b",
    {"albedo_normal": 0.12, "albedo_a": 0.0, "albedo_b": 0.0},
)

#: The envelope, cold to hot.
ENVELOPE: dict[str, AlbedoLaw] = {"cold": BURGER2024, "middle": HAYNE2017, "hot": FOOTE2020}
LAWS: dict[str, AlbedoLaw] = {law.name: law for law in (BURGER2024, HAYNE2017, FOOTE2020, CONSTANT)}

__all__ = ["BURGER2024", "CONSTANT", "ENVELOPE", "FOOTE2020", "HAYNE2017", "LAWS", "AlbedoLaw"]
