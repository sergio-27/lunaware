"""Contract between the three interfaces, and a real agent session over stdio.

The same scenario is run through the Python API, the command line and the MCP server
(started as a subprocess, as an MCP client would). All three must return the same
payload; only wall-clock timings may differ.
"""

from __future__ import annotations

import asyncio
import json
import sys

import pytest
from mcp import Client, StdioServerParameters

from lunarenv import api
from lunarenv.cli import main
from lunarenv.geometry.ephemeris import DE440

pytestmark = pytest.mark.data

REQUIRED = ("pgda78.site04.surf", "pgda90.ldem_80s_80mpp", *DE440.kernel_ids)
SCENARIO = {
    "name": "contract",
    "site": {"x_m": -2822.5, "y_m": -9062.5},
    "period": {
        "start_utc": "2026-06-01T00:00:00",
        "end_utc": "2026-06-08T00:00:00",
        "step_minutes": 120,
    },
    "uncertainty": {"members": 3, "terrain_clones": 1, "seed": 3},
    "illumination": {"settings": {"n_azimuth": 360}},
}


def comparable(payload: dict) -> dict:
    """Payload without wall-clock timings (the only field allowed to differ)."""
    return {k: v for k, v in payload.items() if k != "timing_s"}


def test_api_cli_and_mcp_return_the_same_payload(store, tmp_path, capsys, monkeypatch):
    missing = [i for i in REQUIRED if not store.is_available(i)]
    if missing:
        pytest.skip(f"not in the data cache: {', '.join(missing)}")
    out = tmp_path / "runs"
    monkeypatch.setenv("LUNARENV_DATA_DIR", str(store.root))

    from_api = api.run_summary(api.run_scenario(SCENARIO, output_dir=out, store=store))

    scenario_file = tmp_path / "scenario.json"
    scenario_file.write_text(json.dumps(SCENARIO))
    assert main(["run", str(scenario_file), "-o", str(out), "--json"]) == 0
    from_cli = json.loads(capsys.readouterr().out)

    async def agent_session() -> dict:
        parameters = StdioServerParameters(
            command=sys.executable,
            args=["-m", "lunarenv.mcp_server"],
            env={"LUNARENV_DATA_DIR": str(store.root), "PATH": "/usr/bin:/bin"},
        )
        async with Client(parameters) as client:
            tools = sorted(t.name for t in (await client.list_tools()).tools)
            status = await client.call_tool("get_validation_status")
            run = await client.call_tool(
                "run_scenario", {"scenario": SCENARIO, "output_dir": str(out)}
            )
            assert not run.is_error, run.content
            return {
                "tools": tools,
                "status": status.structured_content,
                "run": run.structured_content,
            }

    session = asyncio.run(agent_session())
    assert session["tools"] == ["describe_models", "get_validation_status", "run_scenario"]
    assert session["status"] == json.loads(json.dumps(api.get_validation_status()))

    expected = json.loads(json.dumps(comparable(from_api)))
    assert comparable(from_cli) == expected
    assert comparable(session["run"]) == expected
    assert from_api["run_id"] == from_cli["run_id"] == session["run"]["run_id"]
    assert (out / from_api["run_id"] / "results.nc").exists()
    # an agent gets what it needs to report responsibly
    assert (
        session["run"]["validity"]["warnings"]
        and "constant_albedo" in session["run"]["model_form_sensitivity"]
    )
