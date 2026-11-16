"""Single-bounce terrain radiation and Earthshine (thermal.terrain_radiation)."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.data.dem import Dem, GridGeometry
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack
from lunarenv.thermal import terrain_radiation as tr
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.surface import plane_normal, surface_normal, view_factors
from lunarenv.validation import bowl

PROJ = SouthPolarStereographic()


def _flat_stack(half_width: float = 400.0, pixel: float = 5.0) -> TerrainStack:
    n = round(2 * half_width / pixel)
    grid = GridGeometry(-half_width, -60000.0 + half_width, pixel, n, n, PROJ)
    dem = Dem(np.zeros((n, n), dtype=np.float32), grid)
    return TerrainStack((TerrainLevel(dem, 0.0, half_width * 1.4, pixel / 2),))


@pytest.fixture(scope="module")
def bowl_stack() -> TerrainStack:
    return bowl.bowl_stack()


def test_flat_ground_level_receiver_sees_no_terrain() -> None:
    view = tr.visible_terrain(_flat_stack(), 0.0, -60000.0, n_rays=72)
    assert tr.cosine_weights(view, [0.0, 0.0, 1.0]).sum() < 1e-6


def test_vertical_wall_on_flat_ground_sees_half_terrain() -> None:
    view = tr.visible_terrain(_flat_stack(), 0.0, -60000.0, height_above_surface_m=2.0, n_rays=144)
    for azimuth in (0.0, 90.0, 225.0):
        weight = tr.cosine_weights(view, plane_normal(np.pi / 2, np.radians(azimuth))).sum()
        # the ground beyond the DEM edge (above -0.3 degrees) is missing
        assert weight == pytest.approx(0.5, abs=5e-3)


@pytest.mark.parametrize("offset", [0.0, 300.0])
def test_terrain_weight_equals_ground_view_factor(bowl_stack: TerrainStack, offset: float) -> None:
    x, y = bowl.CENTRE_XY_M[0] + offset, bowl.CENTRE_XY_M[1]
    view = tr.visible_terrain(bowl_stack, x, y, n_rays=360)
    normal = surface_normal(bowl_stack.levels[0].dem, x, y)
    horizon = compute_horizons(bowl_stack, x, y, n_azimuth=360)
    factors = view_factors(horizon.elevation_rad[0].astype(float), horizon.azimuth_rad, normal)
    assert tr.cosine_weights(view, normal).sum() == pytest.approx(factors.ground, rel=1e-3)


def test_bowl_view_factor_matches_sphere(bowl_stack: TerrainStack) -> None:
    # from the floor centre the cap occupies (2 pi R H) / (4 pi R^2) = H / 2R of the view
    view = tr.visible_terrain(bowl_stack, *bowl.CENTRE_XY_M, n_rays=360)
    expected = bowl.DEPTH_M / (2 * bowl.sphere_radius_m())
    assert tr.cosine_weights(view, [0.0, 0.0, 1.0]).sum() == pytest.approx(expected, rel=5e-3)


def test_patches_need_a_common_point(bowl_stack: TerrainStack) -> None:
    x, y = bowl.CENTRE_XY_M
    a = tr.visible_terrain(bowl_stack, x, y, n_rays=36)
    b = tr.visible_terrain(bowl_stack, x + 10.0, y, n_rays=36)
    with pytest.raises(ValueError, match="share"):
        tr.build_patches([a, b])


def test_patch_points_stay_in_their_patch(bowl_stack: TerrainStack) -> None:
    view = tr.visible_terrain(bowl_stack, *bowl.CENTRE_XY_M, n_rays=72)
    patches, (index,) = tr.build_patches([view], points_per_patch=3)
    assert np.all(np.diff(patches.owner) >= 0)
    assert set(np.unique(patches.owner)) == set(range(patches.size))
    assert np.bincount(patches.owner).max() <= 3
    # every representative point is a visible sample of its own patch
    for i in range(patches.n_points):
        mine = index == patches.owner[i]
        d = np.hypot(view.sample_x_m[mine] - patches.x_m[i], view.sample_y_m[mine] - patches.y_m[i])
        assert d.min() < 1e-6


def test_bowl_default_layout_close_to_analytic() -> None:
    sectors, edges = tr.patch_layout()
    result = bowl.run_bowl(n_sectors=sectors, band_edges_m=edges, n_rays=360, patch_azimuths=360)
    assert all(result.in_shadow)
    for vis, ir in zip(result.visible_w_m2, result.infrared_w_m2, strict=True):
        assert vis == pytest.approx(result.reference_visible_w_m2, rel=0.05)
        assert ir == pytest.approx(result.reference_infrared_w_m2, rel=0.05)


def test_diffuse_albedo_of_constant_law_is_the_constant() -> None:
    props = RegolithProperties(albedo_normal=0.2, albedo_a=0.0, albedo_b=0.0)
    assert tr.diffuse_albedo(props) == pytest.approx(0.2, rel=1e-6)
    # the incidence-dependent law is brighter under diffuse light than at normal incidence
    assert tr.diffuse_albedo(RegolithProperties()) > RegolithProperties().albedo_normal


def test_patch_layout_refinement_doubles_sectors_and_bands() -> None:
    s1, e1 = tr.patch_layout(1)
    s2, e2 = tr.patch_layout(2)
    assert s2 == 2 * s1
    assert len(e2) - 3 == 2 * (len(e1) - 3)  # bands between 25 m and the last finite edge
    assert e1[-2] == pytest.approx(e2[-2])


def test_earthshine_full_and_new_earth() -> None:
    stack = _flat_stack()
    horizons = compute_horizons(stack, 0.0, -60000.0, n_azimuth=72)
    site = PROJ.to_cartesian(0.0, -60000.0, 0.0)
    _, e_y, up = PROJ.local_basis(0.0, -60000.0)
    d = tr.EARTH_MOON_MEAN_DISTANCE_M
    # Earth at 30 degrees elevation; Sun behind the observer (full Earth) or behind the
    # Earth (new Earth)
    direction = np.cos(np.radians(30)) * e_y + np.sin(np.radians(30)) * up
    earth = np.tile(site + d * direction, (2, 1))
    au = 1.496e11
    sun_full = np.tile(site - au * direction, (2, 1))
    sun_new = np.tile(site + au * direction, (2, 1))
    et = np.array([0.0, 3600.0])
    full = tr.earthshine(horizons, et, earth, sun_full, 6.371e6, projection=PROJ)
    new = tr.earthshine(horizons, et, earth, sun_new, 6.371e6, projection=PROJ)
    assert full.illuminated_fraction[0] == pytest.approx(1.0, abs=1e-6)
    assert new.illuminated_fraction[0] == pytest.approx(0.0, abs=1e-6)
    assert full.reflected_normal_w_m2[0] == pytest.approx(0.100, rel=1e-3)
    assert full.infrared_normal_w_m2[0] == pytest.approx(0.056, rel=1e-3)
    assert new.reflected_normal_w_m2[0] == pytest.approx(0.0, abs=1e-6)
    vis, ir, _ = full.on([0.0, 0.0, 1.0])
    assert vis[0] == pytest.approx(0.100 * 0.5, rel=1e-3)  # sin(30 deg)
    assert ir[0] == pytest.approx(0.056 * 0.5, rel=1e-3)
    # a wall facing away from the Earth receives nothing
    vis, ir, _ = full.on(plane_normal(np.pi / 2, np.pi))
    assert vis[0] == 0.0 and ir[0] == 0.0


def test_selected_patches_only() -> None:
    stack = bowl.bowl_stack()
    view = tr.visible_terrain(stack, *bowl.CENTRE_XY_M, n_rays=72)
    patches, _ = tr.build_patches([view])
    light = tr.patch_sunlight(
        stack, patches, np.arange(3.0) * 3600, bowl.sun_positions(), bowl.SUN_RADIUS_M
    )
    only = np.zeros(patches.size, dtype=bool)
    only[: patches.size // 2] = True
    props = RegolithProperties(albedo_a=0.0, albedo_b=0.0)
    some = tr.patch_exitance(
        light, props, step_s=3600.0, spinup_samples=1, equilibrium=True, only=only
    )
    every = tr.patch_exitance(light, props, step_s=3600.0, spinup_samples=1, equilibrium=True)
    assert np.all(some.reflected_w_m2[~only] == 0.0)
    np.testing.assert_allclose(some.reflected_w_m2[only], every.reflected_w_m2[only])
