# MCP server

`lunarenv-mcp` serves the three lunarenv operations to AI agents over the Model Context
Protocol (standard input/output transport). It contains no logic: each tool calls the
function of the same name in `lunarenv.api`, and a test checks that the Python API, the
command line and the server return the same payload for the same scenario.

| Kind | Name | What it does |
|---|---|---|
| tool | `run_scenario(scenario, output_dir="runs")` | Runs a scenario, writes `<output_dir>/<run_id>/` and returns the run summary. |
| tool | `describe_models()` | Assumptions, simplifications, validity range, parameters with uncertainty, numerical error and validation cases of each model. |
| tool | `get_validation_status()` | Stored outcome of every validation case; computes nothing. |
| resource | `lunarenv://scenario-schema` | JSON Schema of a scenario, so that an agent can write a valid one. |

## Client configuration

For Claude Code, in the project's `.mcp.json` (adapt the paths):

```json
{
  "mcpServers": {
    "lunarenv": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/lunarenv", "lunarenv-mcp"],
      "env": { "LUNARENV_DATA_DIR": "/path/to/lunarenv/.data" }
    }
  }
}
```

Any MCP client that can launch a subprocess works the same way. The data cache must
already contain the terrain and kernels (`lunarenv data fetch kernels site04 regional`);
the server never downloads on its own.

## What the agent receives

`run_scenario` returns a summary, not the time series: `run_id`, result checksum, site,
nominal values (`summary`), ensemble quantiles (`uncertainty`), `numerical_error`,
`model_form_sensitivity`, `validity.warnings`, the list of variables with units, timings
and the paths of `results.nc`, `metadata.json` and `scenario.json`.

The server's instructions tell the agent to report uncertainty, validity warnings and
model-form sensitivity together with nominal values, and to consult the validation
status before relying on a result. That is guidance to the client, not something the
server can enforce.

## Behaviour to know about

- **Duration.** A one-year hourly run with default uncertainty takes one to two minutes.
  The computation runs in a worker thread, so the server stays responsive, but the
  client's tool-call timeout must allow for it.
- **Errors.** An invalid scenario, a site outside the terrain tile or missing data come
  back as tool errors with the original message. Anything else is a defect and is
  reported as a generic failure, with the traceback in the server log.
- **Files.** Results are written on the machine that runs the server, under
  `output_dir` relative to the server's working directory unless an absolute path is
  given.
- **One run at a time.** SPICE keeps process-wide state; concurrent runs in one server
  process are refused rather than interleaved.
- **Transport.** Only standard input/output is configured. The SDK's HTTP transports
  are not enabled; exposing the server on a network would need authentication, which
  is out of scope for the MVP.
