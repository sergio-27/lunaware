"""MCP server: it exposes exactly the three API operations and adds nothing of its own."""

from __future__ import annotations

import asyncio
import json

from mcp import Client

from lunarenv import api
from lunarenv.mcp_server import server


def call(coroutine_factory):
    async def session():
        async with Client(server) as client:
            return await coroutine_factory(client)

    return asyncio.run(session())


def test_tools_are_the_three_operations_and_nothing_else():
    tools = call(lambda client: client.list_tools()).tools
    assert sorted(t.name for t in tools) == [
        "describe_models",
        "get_validation_status",
        "run_scenario",
    ]
    by_name = {t.name: t for t in tools}
    assert all(t.description for t in tools)
    schema = by_name["run_scenario"].input_schema
    assert schema["required"] == ["scenario"] and set(schema["properties"]) == {
        "scenario",
        "output_dir",
    }


def test_describe_models_and_validation_status_equal_the_api():
    described = call(lambda client: client.call_tool("describe_models"))
    assert not described.is_error
    assert described.structured_content == {"result": api.describe_models()}
    status = call(lambda client: client.call_tool("get_validation_status"))
    assert status.structured_content == json.loads(json.dumps(api.get_validation_status()))


def test_scenario_schema_resource_equals_the_api():
    resources = call(lambda client: client.list_resources()).resources
    assert [str(r.uri) for r in resources] == ["lunarenv://scenario-schema"]
    content = call(lambda client: client.read_resource("lunarenv://scenario-schema")).contents[0]
    assert json.loads(content.text) == api.json_schema()


def test_an_invalid_scenario_is_reported_as_a_tool_error_not_a_crash():
    result = call(
        lambda client: client.call_tool(
            "run_scenario", {"scenario": {"name": "x", "site": {"x_m": 1.0}}}
        )
    )
    assert result.is_error
    assert "validation error" in result.content[0].text.lower()
    # the server is still alive afterwards
    assert not call(lambda client: client.call_tool("get_validation_status")).is_error


def test_instructions_tell_an_agent_how_to_report():
    assert server.instructions and "validity.warnings" in server.instructions
