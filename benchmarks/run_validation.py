"""Run validation cases, store the results with the package and regenerate VALIDATION.md.

Usage:  LUNARENV_DATA_DIR=.data uv run python benchmarks/run_validation.py [--slow] [CASE_ID ...]

Without ids every fast case runs; ``--slow`` adds the cases that take minutes. Stored
results of cases that are not re-run are kept. (The ``lunarenv validate`` command of
Phase 5 wraps the same calls.)
"""

from __future__ import annotations

import sys
from pathlib import Path

from lunarenv.data.store import DataStore
from lunarenv.validation import load_report, render_markdown, run_cases, save_report

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    arguments = sys.argv[1:]
    slow = "--slow" in arguments
    ids = [a for a in arguments if not a.startswith("--")] or None
    report = run_cases(
        DataStore(), ids, include_slow=slow or ids is not None, previous=load_report()
    )
    save_report(report)
    (ROOT / "VALIDATION.md").write_text(render_markdown(report), encoding="utf-8")
    for result in report.results:
        took = f" ({result.seconds:.0f} s)" if result.seconds else ""
        print(f"{result.id:8s} {result.status!s:8s} {result.summary}{took}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
