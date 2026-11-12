"""MCP server: the three lunarenv operations for AI agents.

The server has no logic of its own. Each tool calls the function of the same name in
:mod:`lunarenv.api` and returns what it returns; the scenario schema is offered as a
resource so that an agent can write a valid scenario. The physics is computed by the
same deterministic code as for a human user: the model is never asked to compute.

Run it over standard input/output with ``lunarenv-mcp`` (see ``docs/mcp.md`` for client
configuration).
"""

from __future__ import annotations

import json
from typing import Any

import anyio
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from lunarenv import __version__, api

INSTRUCTIONS = """\
lunarenv computes the illumination and regolith temperature at a point near the lunar south pole \
over a period of time, with uncertainty and provenance. All numbers come from deterministic code.

Typical use:
1. Read the resource lunarenv://scenario-schema to see every field of a scenario and its units.
2. Call run_scenario with a scenario object. Minimal example:
   {"name": "rim-2026", "site": {"terrain": "site04", "x_m": -2822.5, "y_m": -9062.5},
    "period": {"start_utc": "2026-01-01T00:00:00", "end_utc": "2027-01-01T00:00:00"}}
   A one-year hourly run with the default uncertainty settings takes one to two minutes.
3. Read the returned summary. Always report the uncertainty quantiles, the entries of \
validity.warnings and model_form_sensitivity together with the nominal values: at the pole the \
thermal results are outside the validated range of the regolith model.
4. Call get_validation_status before relying on a result, and describe_models for the assumptions.

Arrays are written to results.nc (NetCDF) and full provenance to metadata.json in the returned \
directory; the tool result is a summary, not the time series.
"""

server = MCPServer("lunarenv", instructions=INSTRUCTIONS, version=__version__)


@server.tool(
    description=(
        "Run a lunarenv scenario (illumination, regolith temperature, uncertainty) and write the "
        "results to <output_dir>/<run_id>/. Returns a summary with run_id, result checksum, "
        "nominal values, uncertainty quantiles, numerical error, model-form sensitivity, validity "
        "warnings and the paths of the result files. The same scenario and data always give the "
        "same run_id."
    )
)
async def run_scenario(scenario: dict[str, Any], output_dir: str = "runs") -> dict[str, Any]:
    def work() -> dict[str, Any]:
        try:
            return api.run_summary(api.run_scenario(scenario, output_dir=output_dir))
        except (ValueError, RuntimeError, FileNotFoundError) as error:
            # invalid scenario, site outside the terrain, data not in the cache: errors the
            # caller can act on are returned as tool errors with the original message
            raise ToolError(str(error)) from error

    # computed in a worker thread so that the server keeps answering while a run is in progress
    return await anyio.to_thread.run_sync(work)


@server.tool(
    description=(
        "Describe every model: assumptions, simplifications, validity range, parameters with their "
        "uncertainty intervals, numerical error estimates, validation cases and references."
    )
)
def describe_models() -> list[dict[str, Any]]:
    return api.describe_models()


@server.tool(
    description=(
        "Stored outcome of every validation case (PASS, FAIL, PARTIAL, PENDING) with its source, "
        "metric, tolerance and result, plus counts by status. Does not recompute anything."
    )
)
def get_validation_status() -> dict[str, Any]:
    return dict(api.get_validation_status())


@server.resource(
    "lunarenv://scenario-schema",
    name="scenario-schema",
    description="JSON Schema of a lunarenv scenario: every field, its unit, default and limits.",
    mime_type="application/json",
)
def scenario_schema() -> str:
    return json.dumps(api.json_schema(), indent=1)


def main() -> None:
    """Serve over standard input/output."""
    server.run("stdio")


if __name__ == "__main__":
    main()
