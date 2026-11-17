"""Command line interface: ``lunarenv run | report | validate | describe | schema | data``.

Thin wrapper over :mod:`lunarenv.api`; it holds no logic of its own. Results go to
standard output as text, or as JSON with ``--json`` (for scripts and agents).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from lunarenv import __version__


def _print_json(value: Any) -> None:
    json.dump(value, sys.stdout, indent=1)
    sys.stdout.write("\n")


def _run(args: argparse.Namespace) -> int:
    from lunarenv.api import run_scenario, run_summary

    result = run_scenario(args.scenario, output_dir=args.output)
    meta = result.metadata
    if args.json:
        _print_json(run_summary(result))
        return 0
    print(f"run {result.run_id}  ->  {result.directory}")
    print(f"result checksum  {result.result_checksum}")
    ill = meta["summary"]["illumination"]
    print(
        f"illumination: mean visible solar disc {100 * ill['mean_visible_fraction']:.2f} %, "
        f"longest shadow {ill['longest_shadow_hours']:.0f} h, longest illumination {ill['longest_illumination_hours']:.0f} h"
    )
    ensemble = meta["uncertainty"].get("illumination")
    if ensemble:
        q = ensemble["mean_visible_fraction"]["quantiles"]
        print(
            f"              terrain ensemble 5/50/95 %: {100 * q[0]:.2f} / {100 * q[1]:.2f} / {100 * q[2]:.2f} %"
        )
    regional = meta["uncertainty"].get("regional_error_bound")
    if regional:
        lo, hi = regional["raised"], regional["lowered"]
        print(
            f"              regional DEM error bound (+/-{regional['k']:g} rms): "
            f"{100 * lo['mean_visible_fraction']:.2f} - {100 * hi['mean_visible_fraction']:.2f} %, "
            f"longest shadow {hi['longest_shadow_hours']:.0f} - {lo['longest_shadow_hours']:.0f} h"
        )
    thermal = meta["summary"].get("thermal")
    if thermal:
        bounds = meta["uncertainty"].get("thermal_envelope", {}).get("bounds")
        if bounds:
            lo_min, hi_min = bounds["surface_min_k"]
            lo_max, hi_max = bounds["surface_max_k"]
            lo_mean, hi_mean = bounds["surface_mean_k"]
            print(
                "thermal bounds (albedo envelope x ensemble 5-95 %): surface minimum "
                f"{lo_min:.1f}-{hi_min:.1f} K, maximum {lo_max:.1f}-{hi_max:.1f} K, "
                f"mean {lo_mean:.1f}-{hi_mean:.1f} K"
            )
        envelope = thermal.get("albedo_envelope")
        if envelope:
            print(
                f"         nominal terrain: cold law max {envelope['cold']['surface_max_k']:.1f} K, "
                f"scenario law {envelope['nominal']['surface_max_k']:.1f} K, "
                f"hot law {envelope['hot']['surface_max_k']:.1f} K"
            )
        print(
            f"thermal, scenario law: surface {thermal['surface_min_k']:.1f} to {thermal['surface_max_k']:.1f} K "
            f"(mean {thermal['surface_mean_k']:.1f} K) on a {thermal['surface_slope_deg']:.2f} deg slope"
        )
        spread = meta["uncertainty"].get("thermal")
        if spread:
            q = spread["surface_max_k"]["quantiles"]
            print(
                f"         surface maximum, ensemble 5/50/95 %: {q[0]:.1f} / {q[1]:.1f} / {q[2]:.1f} K"
            )
        alternative = meta["model_form_sensitivity"].get("constant_albedo")
        if alternative:
            print(
                f"         with constant albedo: mean {alternative['surface_mean_k']:.1f} K "
                f"({alternative['difference_surface_mean_k']:+.1f} K); not settled by validation"
            )
        k_alt = meta["model_form_sensitivity"].get("conductivity")
        if k_alt:
            changes = ", ".join(
                f"{d} {v:+.1f} K" for d, v in k_alt["difference_depth_mean_k"].items()
            )
            print(
                f"         with {k_alt['conductivity_model']} conductivity: minimum "
                f"{k_alt['difference_surface_min_k']:+.1f} K, mean "
                f"{k_alt['difference_surface_mean_k']:+.1f} K; means at depth {changes}"
            )
        grid = thermal.get("slope_grid")
        if grid:
            lo, hi = grid["range_surface_max_k"]
            slopes = sorted({o["slope_deg"] for o in grid["orientations"]})
            print(
                f"         slope grid {slopes[0]:g}-{slopes[-1]:g} deg, all aspects and laws: "
                f"surface maximum {lo:.1f} to {hi:.1f} K"
            )
        diffuse = thermal.get("diffuse")
        if diffuse:
            nominal_case = diffuse["cases"]["nominal"]
            print(
                f"         diffuse at the minimum: terrain {nominal_case['terrain_irradiance_at_minimum_w_m2']:.2g} "
                f"W m-2, Earthshine absorbed {nominal_case['earthshine_absorbed_at_minimum_w_m2']:.2g} W m-2; "
                f"minimum raised {nominal_case['terrain_effect_k']['surface_min']:+.2f} K (terrain), "
                f"{nominal_case['earthshine_effect_k']['surface_min']:+.2f} K (Earthshine)"
            )
    surfaces = meta["summary"].get("oriented_surfaces")
    if surfaces:
        for name, s in surfaces["surfaces"].items():
            print(
                f"surface {name} (tilt {s['tilt_deg']:g}, facing {s['facing_grid_azimuth_deg']:g}): "
                f"direct sun mean {s['mean_direct_irradiance_w_m2']:.0f} W m-2, max "
                f"{s['max_direct_irradiance_w_m2']:.0f} W m-2, lit {100 * s['fraction_of_time_with_direct_sun']:.0f} %; "
                f"view factor sky {s['sky_view_factor']:.3f}, ground {s['ground_view_factor']:.3f}"
            )
            if "mean_terrain_scattered_w_m2" in s:
                print(
                    f"         diffuse mean: terrain {s['mean_terrain_scattered_w_m2']:.2g} W m-2 scattered, "
                    f"{s['mean_terrain_infrared_w_m2']:.2g} W m-2 infrared; Earthshine "
                    f"{s['mean_earthshine_reflected_w_m2'] + s['mean_earthshine_infrared_w_m2']:.2g} W m-2"
                )
    for warning in meta["validity"]["warnings"]:
        print(f"WARNING {warning}")
    print(f"wall time {meta['timing_s']['total_s']:.1f} s")
    return 0


def _validate(args: argparse.Namespace) -> int:
    import textwrap

    from lunarenv.validation import core

    not_rerun: dict[str, str] = {}
    if args.run:
        from lunarenv.data.store import DataStore

        store = DataStore()
        print(f"data cache: {store.root}", file=sys.stderr)
        report = core.run_cases(
            store,
            args.cases or None,
            include_slow=args.slow or bool(args.cases),
            previous=core.load_report(),
        )
        core.save_report(report)
        not_rerun = {k: v for k, v in report.not_rerun.items() if v != "not selected"}
        if args.markdown:
            Path(args.markdown).write_text(core.render_markdown(report), encoding="utf-8")
    status = core.get_validation_status()
    if args.json:
        _print_json(status)
    else:
        print(
            f"lunarenv {status['lunarenv_version']}, last validation run {status['generated_utc']}"
        )
        for case in status["cases"]:
            computed = (case["run_utc"] or "-")[:10]
            summary = textwrap.shorten(case["summary"], width=96, placeholder=" ...")
            print(f"{case['id']:8s} {case['status']:8s} {computed:10s} {summary}")
        print("totals: " + ", ".join(f"{n} {s}" for s, n in sorted(status["counts"].items())))
    if not_rerun:
        sys.stdout.flush()
        print(
            "\nnot re-run in this run (any stored result above is from an earlier run):",
            file=sys.stderr,
        )
        for case_id, reason in not_rerun.items():
            print(f"  {case_id:8s} {reason}", file=sys.stderr)
    return 1 if "FAIL" in status["counts"] and args.strict else 0


def _describe(args: argparse.Namespace) -> int:
    from lunarenv.api import describe_models

    models = describe_models()
    if args.json:
        _print_json(models)
        return 0
    for model in models:
        print(f"{model['name']} {model['version']}\n  {model['summary']}")
        for heading in ("assumptions", "simplifications", "validity_range", "uncertainty_sources"):
            print(f"  {heading.replace('_', ' ')}:")
            for item in model[heading]:
                print(f"    - {item}")
        print(f"  numerical error: {model['numerical_error']}")
        print(f"  validation cases: {', '.join(model['validation_cases'])}\n")
    return 0


def _report(args: argparse.Namespace) -> int:
    from lunarenv.report import write_report

    run = args.run if args.run.is_dir() else Path("runs") / args.run
    print(write_report(run, fmt=args.format, output=args.output))
    return 0


def _export(args: argparse.Namespace) -> int:
    from lunarenv.export import export_run

    run = args.run if args.run.is_dir() else Path("runs") / args.run
    manifest = export_run(
        run,
        args.output,
        mesh_half_width_m=args.mesh_half_width,
        mesh_stride=args.mesh_stride,
        ring_radius_m=args.ring_radius,
    )
    if args.json:
        _print_json(manifest)
        return 0
    print(manifest["directory"])
    for entry in manifest["files"]:
        size = f"{entry['rows']} rows" if "rows" in entry else f"{entry['faces']} faces"
        print(f"  {entry['name']}: {size}, {entry['bytes']} bytes")
    print("  manifest.json: frame, units, checksums, validation status")
    return 0


def _sites(args: argparse.Namespace) -> int:
    from lunarenv.sites import site_status

    rows = site_status()
    if args.json:
        _print_json(rows)
        return 0
    print(f"{'terrain':<10s} {'DEM':<12s} {'clones':<7s} description")
    for row in rows:
        state = "cached" if row["cached"] else "pinned" if row["pinned"] else "not pinned"
        clones = "-" if row["clones_cached"] is None else str(row["clones_cached"])
        print(f"{row['name']:<10s} {state:<12s} {clones:<7s} {row['description']}")
    print(
        "cached: ready to run. pinned: `lunarenv data fetch <terrain>`. "
        "not pinned: a maintainer fetches it once with --pin."
    )
    return 0


def _schema(_: argparse.Namespace) -> int:
    from lunarenv.api import json_schema

    _print_json(json_schema())
    return 0


def _data(args: argparse.Namespace) -> int:
    from lunarenv.data.__main__ import main as data_main

    return data_main(args.arguments)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lunarenv", description="Open, verified lunar surface environment engine."
    )
    parser.add_argument("--version", action="version", version=f"lunarenv {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run a scenario (YAML or JSON)")
    run.add_argument("scenario", type=Path)
    run.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("runs"),
        help="directory that receives <run_id>/ (default: runs)",
    )
    run.add_argument("--json", action="store_true", help="print a JSON summary")
    run.set_defaults(handler=_run)

    validate = sub.add_parser("validate", help="show the stored validation status, or re-run cases")
    validate.add_argument("cases", nargs="*", help="case ids to re-run (with --run)")
    validate.add_argument(
        "--run", action="store_true", help="re-run cases and update the stored results"
    )
    validate.add_argument(
        "--slow", action="store_true", help="with --run: include cases that take minutes"
    )
    validate.add_argument(
        "--markdown", metavar="FILE", help="with --run: also write the report as Markdown"
    )
    validate.add_argument(
        "--strict", action="store_true", help="exit with status 1 if any case is FAIL"
    )
    validate.add_argument("--json", action="store_true")
    validate.set_defaults(handler=_validate)

    describe = sub.add_parser(
        "describe", help="assumptions, validity and uncertainty of each model"
    )
    describe.add_argument("--json", action="store_true")
    describe.set_defaults(handler=_describe)

    report = sub.add_parser("report", help="write the HTML or Markdown report of a run")
    report.add_argument("run", type=Path, help="run directory, or a run id under ./runs")
    report.add_argument("--format", choices=("html", "md"), default="html")
    report.add_argument(
        "-o", "--output", type=Path, help="output file (default: inside the run directory)"
    )
    report.set_defaults(handler=_report)

    export = sub.add_parser(
        "export",
        help="write a run as CSV tables and OBJ meshes for another thermal tool",
    )
    export.add_argument("run", type=Path, help="run directory, or a run id under ./runs")
    export.add_argument(
        "-o", "--output", type=Path, help="output directory (default: <run>/export)"
    )
    export.add_argument(
        "--ring-radius",
        type=float,
        metavar="M",
        help="also write horizon_ring.obj: the horizon as a wall of this radius (m)",
    )
    export.add_argument(
        "--mesh-half-width",
        type=float,
        metavar="M",
        help="also write terrain.obj: the site DEM within this distance (m) in x and y",
    )
    export.add_argument(
        "--mesh-stride",
        type=int,
        default=1,
        metavar="N",
        help="with --mesh-half-width: keep every N-th DEM pixel (default 1)",
    )
    export.add_argument("--json", action="store_true", help="print the manifest")
    export.set_defaults(handler=_export)

    sites = sub.add_parser("sites", help="list the terrain sets and which can be run")
    sites.add_argument("--json", action="store_true")
    sites.set_defaults(handler=_sites)

    schema = sub.add_parser("schema", help="print the JSON Schema of a scenario")
    schema.set_defaults(handler=_schema)

    data = sub.add_parser("data", help="manage the data cache: fetch | verify | list")
    data.add_argument("arguments", nargs=argparse.REMAINDER)
    data.set_defaults(handler=_data)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.handler(args))
    except (ValueError, FileNotFoundError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
