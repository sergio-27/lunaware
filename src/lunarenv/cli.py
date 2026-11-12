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
    thermal = meta["summary"].get("thermal")
    if thermal:
        print(
            f"thermal: surface {thermal['surface_min_k']:.1f} to {thermal['surface_max_k']:.1f} K "
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
    for warning in meta["validity"]["warnings"]:
        print(f"WARNING {warning}")
    print(f"wall time {meta['timing_s']['total_s']:.1f} s")
    return 0


def _validate(args: argparse.Namespace) -> int:
    from lunarenv.validation import core

    if args.run:
        from lunarenv.data.store import DataStore

        report = core.run_cases(
            DataStore(),
            args.cases or None,
            include_slow=args.slow or bool(args.cases),
            previous=core.load_report(),
        )
        core.save_report(report)
        if args.markdown:
            Path(args.markdown).write_text(core.render_markdown(report), encoding="utf-8")
    status = core.get_validation_status()
    if args.json:
        _print_json(status)
    else:
        print(
            f"lunarenv {status['lunarenv_version']}, validation results of {status['generated_utc']}"
        )
        for case in status["cases"]:
            print(f"{case['id']:8s} {case['status']:8s} {case['summary'][:110]}")
        print("totals: " + ", ".join(f"{n} {s}" for s, n in sorted(status["counts"].items())))
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
