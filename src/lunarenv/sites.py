"""Terrain sets: which catalog products describe the surroundings of a site."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TerrainSet:
    """Catalog ids of the terrain around a site."""

    name: str
    description: str
    near_id: str
    regional_id: str
    clone_group: str
    regional_error_id: str = "pgda90.ldem_80s_80mpp_err"

    @property
    def required_ids(self) -> tuple[str, ...]:
        return (self.near_id, self.regional_id)


TERRAIN_SETS: dict[str, TerrainSet] = {
    "site04": TerrainSet(
        "site04",
        "Shackleton rim (PGDA Site 4): 5 m site DEM, 80 m regional DEM poleward of 80 S",
        "pgda78.site04.surf",
        "pgda90.ldem_80s_80mpp",
        "site04-clones",
    ),
    "site01": TerrainSet(
        "site01",
        "Connecting ridge (PGDA Site 1): 5 m site DEM, 80 m regional DEM poleward of 80 S",
        "pgda78.site01.surf",
        "pgda90.ldem_80s_80mpp",
        "site01-clones",
    ),
}
