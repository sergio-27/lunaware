"""Validation framework.

A *case* compares something the code computes with a reference that the code did not
produce: a closed-form solution, an independent implementation, or published data. Each
case states its source, its metric and its tolerance; tolerances are fixed in the case
definition, before the comparison is run, and changing one is a reviewed code change.

Statuses
    PASS      every check with a tolerance is inside it
    FAIL      at least one check is outside its tolerance (reported, never hidden)
    PARTIAL   the comparison ran on incomplete inputs or covers only part of the claim
    PENDING   the reference could not be obtained; the case says what is missing
    INFO      reported for context, no pass/fail claim

Results are stored with the package (``results.json``) together with the code version
and the checksums of the data they used, because some cases take minutes. That file is
what ``get_validation_status`` reports and what ``VALIDATION.md`` is generated from; it
is regenerated with ``lunarenv validate`` (Phase 5) or ``run_cases``.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from importlib import resources
from pathlib import Path
from typing import Literal, TypedDict

from pydantic import BaseModel, ConfigDict, Field

from lunarenv import __version__
from lunarenv.data.store import DataStore

_RESULTS_FILE = "results.json"


class Status(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    PARTIAL = "PARTIAL"
    PENDING = "PENDING"
    INFO = "INFO"
    NOT_RUN = "NOT RUN"


class Check(BaseModel):
    """One compared quantity."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    quantity: str
    reference: float | None = None
    measured: float | None = None
    tolerance: float | None = Field(None, description="Allowed |measured - reference|.")
    unit: str = ""
    informational: bool = False

    @property
    def passed(self) -> bool | None:
        if self.informational or self.tolerance is None:
            return None
        if self.reference is None or self.measured is None:
            return None
        return abs(self.measured - self.reference) <= self.tolerance


class CaseResult(BaseModel):
    """Outcome of one validation case."""

    model_config = ConfigDict(extra="forbid")

    id: str
    title: str
    kind: Literal["analytic", "property", "cross-check", "published"]
    source: str
    metric: str
    tolerance: str
    status: Status
    summary: str
    checks: list[Check] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)
    data: dict[str, str] = Field(default_factory=dict, description="catalog id -> sha256 used")
    seconds: float | None = None
    run_utc: str | None = Field(
        None, description="When this case was computed; None if not recorded (older results)."
    )


class ValidationReport(BaseModel):
    """Stored results of a validation run."""

    model_config = ConfigDict(extra="forbid")

    schema_version: int = 1
    lunarenv_version: str = __version__
    generated_utc: str = ""
    results: list[CaseResult] = Field(default_factory=list)
    not_rerun: dict[str, str] = Field(
        default_factory=dict,
        description="Cases this run did not compute, with the reason; their entry, if any, "
        "is the stored result of an earlier run (see its run_utc).",
    )

    def by_id(self) -> dict[str, CaseResult]:
        return {r.id: r for r in self.results}


@dataclass(frozen=True, slots=True)
class Case:
    """Definition of a validation case.

    ``run`` receives the data store and returns the checks, a summary, notes and an
    optional status override (``None`` means: PASS if all toleranced checks pass, else
    FAIL). A case without ``run`` is PENDING and must say why.
    """

    id: str
    title: str
    kind: Literal["analytic", "property", "cross-check", "published"]
    source: str
    metric: str
    tolerance: str
    requires: tuple[str, ...] = ()
    slow: bool = False
    run: Callable[[DataStore], Outcome] | None = None
    pending_reason: str = ""


@dataclass(frozen=True, slots=True)
class Outcome:
    checks: list[Check]
    summary: str
    notes: tuple[str, ...] = ()
    status: Status | None = None


def _status_from_checks(checks: Iterable[Check]) -> Status:
    verdicts = [c.passed for c in checks if c.passed is not None]
    if not verdicts:
        return Status.INFO
    return Status.PASS if all(verdicts) else Status.FAIL


def registry() -> tuple[Case, ...]:
    """All cases, in report order."""
    from lunarenv.validation import cases_illumination, cases_thermal

    return (*cases_illumination.CASES, *cases_thermal.CASES)


def run_cases(
    store: DataStore,
    ids: Iterable[str] | None = None,
    *,
    include_slow: bool = True,
    previous: ValidationReport | None = None,
) -> ValidationReport:
    """Run validation cases and return a report.

    Cases that are skipped (not selected, slow and excluded, or missing data) keep their
    entry from ``previous`` when there is one, so a quick run does not erase the stored
    result of an expensive case.
    """
    import time

    selected = set(ids) if ids is not None else None
    known = previous.by_id() if previous else {}
    results: list[CaseResult] = []
    not_rerun: dict[str, str] = {}
    now = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    for case in registry():
        base = {
            "id": case.id,
            "title": case.title,
            "kind": case.kind,
            "source": case.source,
            "metric": case.metric,
            "tolerance": case.tolerance,
        }
        if case.run is None:
            results.append(CaseResult(**base, status=Status.PENDING, summary=case.pending_reason))
            continue
        skip = (selected is not None and case.id not in selected) or (
            case.slow and not include_slow
        )
        missing = [r for r in case.requires if not store.is_available(r)]
        if skip or missing:
            if selected is not None and case.id not in selected:
                reason = "not selected"
            elif case.slow and not include_slow:
                reason = "slow case; include it with --slow"
            else:
                reason = f"data not in the cache {store.root}: {', '.join(missing)}"
            not_rerun[case.id] = reason
            if case.id in known and known[case.id].status is not Status.NOT_RUN:
                results.append(known[case.id])
            else:
                results.append(CaseResult(**base, status=Status.NOT_RUN, summary=reason))
            continue
        start = time.perf_counter()
        outcome = case.run(store)
        results.append(
            CaseResult(
                **base,
                status=outcome.status or _status_from_checks(outcome.checks),
                summary=outcome.summary,
                checks=outcome.checks,
                notes=list(outcome.notes),
                data={r: store.record(r).sha256 for r in case.requires},
                seconds=time.perf_counter() - start,
                run_utc=now,
            )
        )
    return ValidationReport(generated_utc=now, results=results, not_rerun=not_rerun)


def results_path() -> Path:
    return Path(str(resources.files("lunarenv.validation").joinpath(_RESULTS_FILE)))


def load_report(path: Path | None = None) -> ValidationReport | None:
    path = path or results_path()
    if not path.exists():
        return None
    return ValidationReport.model_validate_json(path.read_text(encoding="utf-8"))


def save_report(report: ValidationReport, path: Path | None = None) -> None:
    path = path or results_path()
    path.write_text(json.dumps(report.model_dump(mode="json"), indent=1) + "\n", encoding="utf-8")


class CaseStatus(TypedDict):
    id: str
    title: str
    kind: str
    source: str
    metric: str
    tolerance: str
    status: str
    summary: str
    run_utc: str | None


class ValidationStatus(TypedDict):
    lunarenv_version: str
    generated_utc: str | None
    counts: dict[str, int]
    cases: list[CaseStatus]


def get_validation_status(path: Path | None = None) -> ValidationStatus:
    """Stored validation status: one entry per case plus counts by status.

    Reads the results shipped with the package; it does not recompute anything.
    """
    report = load_report(path)
    stored = report.by_id() if report else {}
    cases: list[CaseStatus] = []
    for case in registry():
        result = stored.get(case.id)
        cases.append(
            CaseStatus(
                id=case.id,
                title=case.title,
                kind=case.kind,
                source=case.source,
                metric=case.metric,
                tolerance=case.tolerance,
                status=str(result.status if result else Status.NOT_RUN),
                summary=result.summary if result else "no stored result",
                run_utc=result.run_utc if result else None,
            )
        )
    counts: dict[str, int] = {}
    for entry in cases:
        counts[entry["status"]] = counts.get(entry["status"], 0) + 1
    return ValidationStatus(
        lunarenv_version=report.lunarenv_version if report else __version__,
        generated_utc=report.generated_utc if report else None,
        counts=counts,
        cases=cases,
    )


def _format(value: float | None, unit: str) -> str:
    if value is None:
        return "-"
    text = f"{value:.4g}" if abs(value) < 1e4 else f"{value:.6g}"
    return f"{text} {unit}".strip()


def _cell(text: str) -> str:
    """Text made safe for a Markdown table cell."""
    return text.replace("|", "\\|").replace("\n", " ")


def render_markdown(report: ValidationReport) -> str:
    """VALIDATION.md: summary table, then the checks and notes of every case."""
    lines = [
        "# Validation status",
        "",
        "<!-- generated by lunarenv.validation.render_markdown; do not edit by hand -->",
        "",
        f"lunarenv {report.lunarenv_version}; results generated {report.generated_utc}.",
        "Reference keys resolve in [docs/references.md](docs/references.md). Tolerances are part of",
        "each case definition (`src/lunarenv/validation/`). They were fixed before the comparison ran;",
        "where a tolerance was written down after exploratory runs, the case notes say so.",
        "",
        "Status: **PASS** inside tolerance; **FAIL** outside (reported as found); **PARTIAL** ran on",
        "incomplete inputs or supports only part of the claim; **PENDING** reference not obtained;",
        "**INFO** context without a pass/fail claim.",
        "",
        "| Case | What is compared | Source | Metric | Tolerance | Result | Status |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in report.results:
        lines.append(
            f"| {r.id} | {_cell(r.title)} | {_cell(r.source)} | {_cell(r.metric)} | {_cell(r.tolerance)} "
            f"| {_cell(r.summary)} | **{r.status}** |"
        )
    counts: dict[str, int] = {}
    for r in report.results:
        counts[str(r.status)] = counts.get(str(r.status), 0) + 1
    lines += ["", "Totals: " + ", ".join(f"{n} {s}" for s, n in sorted(counts.items())) + ".", ""]
    lines.append("## Details")
    for r in report.results:
        run = f" Computed {r.run_utc}." if r.run_utc else ""
        lines += ["", f"### {r.id} - {r.title}", "", f"Status: **{r.status}**.{run} {r.summary}"]
        if r.checks:
            lines += [
                "",
                "| Quantity | Reference | Measured | Difference | Tolerance | Verdict |",
                "|---|---|---|---|---|---|",
            ]
            for c in r.checks:
                difference = (
                    None if c.reference is None or c.measured is None else c.measured - c.reference
                )
                verdict = {True: "pass", False: "FAIL", None: "info"}[c.passed]
                lines.append(
                    f"| {_cell(c.quantity)} | {_format(c.reference, c.unit)} | {_format(c.measured, c.unit)} "
                    f"| {_format(difference, c.unit)} | {_format(c.tolerance, c.unit)} | {verdict} |"
                )
        for note in r.notes:
            lines += ["", f"- {note}"]
        if r.data:
            lines += [
                "",
                "Data: " + ", ".join(f"`{k}` ({v[:12]})" for k, v in r.data.items()) + ".",
            ]
    return "\n".join(lines) + "\n"
