"""Terrain sets: which catalog products describe the surroundings of a site.

One set per site of PGDA product 78 (5 m LOLA site models of the lunar south pole,
Barker et al. 2021): the site DEM for the near field and the regional 80 m DEM beyond.
A scenario names its set in ``site.terrain``; the observer must lie inside the site DEM.

A set can be run once its site DEM is in the data cache, which needs a pinned checksum
(``lunarenv sites`` shows which are). Terrain error clones are in the catalog for
``site04`` and ``site01`` only; runs at the other sites do not vary the near terrain.
"""

from __future__ import annotations

from dataclasses import dataclass

_REGIONAL = "pgda90.ldem_80s_80mpp"


@dataclass(frozen=True, slots=True)
class TerrainSet:
    """Catalog ids of the terrain around a site."""

    name: str
    description: str
    near_id: str
    regional_id: str
    clone_group: str | None = None
    regional_error_id: str = "pgda90.ldem_80s_80mpp_err"

    @property
    def required_ids(self) -> tuple[str, ...]:
        return (self.near_id, self.regional_id)


#: PGDA 78 sites: (set name, name in the product, description on the product page).
PGDA78_SITES: tuple[tuple[str, str, str], ...] = (
    ("site01", "Site01", "Connecting ridge"),
    ("site04", "Site04", "Shackleton rim"),
    ("site06", "Site06", "Nobile rim 1"),
    ("site07", "Site07", "Peak near Shackleton"),
    ("site11", "Site11", "de Gerlache rim"),
    ("site20", "Site20", "Leibnitz beta plateau"),
    ("site20v2", "Site20v2", "Leibnitz beta plateau, extended boundaries"),
    ("site23", "Site23", "Malapert massif"),
    ("site42", "Site42", "de Gerlache-Kocher massif"),
    ("haworth", "Haworth", "Haworth"),
    ("shoemaker", "Shoemaker", "Shoemaker"),
    ("dm1", "DM1", "Amundsen rim"),
    ("dm2", "DM2", "Nobile rim 2"),
    ("sl2", "SL2", "de Gerlache rim"),
    ("sl3", "SL3", "Connecting ridge extension"),
    ("npa", "NPA", "Cabeus exterior wall 1"),
    ("npb", "NPB", "Amundsen 1"),
    ("npc", "NPC", "Idel'son L crater 1"),
    ("npd", "NPD", "Malapert crater 1"),
    ("lm1", "LM1", "Shackleton Rim B"),
    ("lm2", "LM2", "Shoemaker Rim A"),
    ("lm3", "LM3", "Shoemaker Rim B"),
    ("lm4", "LM4", "Shoemaker Rim C"),
    ("lm5", "LM5", "Shoemaker Rim D"),
    ("lm6", "LM6", "Shoemaker Rim E"),
    ("lm7", "LM7", "Faustini Rim A"),
    ("lm8", "LM8", "Shoemaker Rim F"),
)

#: Sets whose terrain error clones are in the catalog.
_CLONE_GROUPS = {"site04": "site04-clones", "site01": "site01-clones"}

TERRAIN_SETS: dict[str, TerrainSet] = {
    name: TerrainSet(
        name,
        f"{description} (PGDA {product}): 5 m site DEM, 80 m regional DEM poleward of 80 S",
        f"pgda78.{name}.surf",
        _REGIONAL,
        _CLONE_GROUPS.get(name),
    )
    for name, product, description in PGDA78_SITES
}

DEFAULT_TERRAIN = "site04"


def site_status(store: object | None = None) -> list[dict[str, object]]:
    """One row per terrain set: what it is and whether it can be run from this cache.

    ``pinned``: the site DEM has a pinned checksum, so it can be fetched. ``cached``:
    the site DEM and the regional DEM are in the data cache, so a scenario can run.
    ``clones_cached`` counts the terrain error clones in the cache (``None``: the
    catalog has none for the set).
    """
    from lunarenv.data.catalog import load_pins
    from lunarenv.data.store import DataStore

    data = store if isinstance(store, DataStore) else DataStore()
    pins = load_pins()
    rows: list[dict[str, object]] = []
    for terrain in TERRAIN_SETS.values():
        clones = (
            None
            if terrain.clone_group is None
            else sum(data.is_available(e.id) for e in data.catalog.group(terrain.clone_group))
        )
        rows.append(
            {
                "name": terrain.name,
                "description": terrain.description,
                "near_id": terrain.near_id,
                "pinned": terrain.near_id in pins,
                "cached": all(data.is_available(i) for i in terrain.required_ids),
                "clones_cached": clones,
            }
        )
    return rows
