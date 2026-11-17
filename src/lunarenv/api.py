"""Public Python API.

Three operations, the same ones exposed by the command line and by the MCP server:

* :func:`run_scenario` - run a declarative scenario;
* :func:`describe_models` - what each model assumes, simplifies and has been validated for;
* :func:`get_validation_status` - the stored outcome of every validation case.
"""

from __future__ import annotations

from typing import Any

from lunarenv.export import export_run
from lunarenv.runner import RunResult, describe_models, load_run, run_scenario
from lunarenv.scenario.schema import Scenario, json_schema, load_scenario
from lunarenv.validation.core import get_validation_status

#: Parts of the metadata returned as the summary of a run (see :func:`run_summary`).
_SUMMARY_KEYS = (
    "run_id",
    "result_checksum_sha256",
    "lunarenv_version",
    "site",
    "summary",
    "uncertainty",
    "model_form_sensitivity",
    "numerical_error",
    "validity",
    "variables",
    "timing_s",
)


def run_summary(result: RunResult) -> dict[str, Any]:
    """JSON-serializable summary of a run: what the command line and the MCP server return.

    The full metadata (models, data provenance, environment, validation status) is in
    ``metadata.json`` in the run directory, and the arrays in ``results.nc``.
    """
    summary = {key: result.metadata[key] for key in _SUMMARY_KEYS}
    directory = result.directory
    summary["directory"] = None if directory is None else str(directory)
    summary["files"] = (
        None
        if directory is None
        else {
            "results": str(directory / "results.nc"),
            "metadata": str(directory / "metadata.json"),
            "scenario": str(directory / "scenario.json"),
        }
    )
    return summary


__all__ = [
    "RunResult",
    "Scenario",
    "describe_models",
    "export_run",
    "get_validation_status",
    "json_schema",
    "load_run",
    "load_scenario",
    "run_scenario",
    "run_summary",
]
