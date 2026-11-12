"""Report generation from a real run: content, self-containment, both formats, CLI."""

from __future__ import annotations

import re

import pytest

from lunarenv.api import run_scenario
from lunarenv.cli import main
from lunarenv.geometry.ephemeris import DE440
from lunarenv.report import TOKENS, key_values, write_report

pytestmark = pytest.mark.data

REQUIRED = ("pgda78.site04.surf", "pgda90.ldem_80s_80mpp", *DE440.kernel_ids)
BASE = {
    "name": "report-test",
    "description": "A <b>scenario</b> & its report.",
    "site": {"x_m": -2822.5, "y_m": -9062.5},
    "period": {
        "start_utc": "2026-03-01T00:00:00",
        "end_utc": "2026-03-11T00:00:00",
        "step_minutes": 120,
    },
    "illumination": {
        "settings": {"n_azimuth": 360},
        "map": {"half_width_m": 10.0, "n_azimuth": 90},
    },
}


@pytest.fixture(scope="module")
def runs(store, tmp_path_factory):
    missing = [i for i in REQUIRED if not store.is_available(i)]
    if missing:
        pytest.skip(f"not in the data cache: {', '.join(missing)}")
    out = tmp_path_factory.mktemp("runs")
    full = run_scenario(
        BASE | {"uncertainty": {"members": 3, "terrain_clones": 1}}, output_dir=out, store=store
    )
    bare = run_scenario(
        BASE
        | {
            "name": "bare",
            "thermal": {"enabled": False},
            "uncertainty": {"members": 0, "terrain_clones": 0},
        },
        output_dir=out,
        store=store,
    )
    return full, bare


def test_html_report_is_complete_and_self_contained(runs):
    full, _ = runs
    path = write_report(full.directory)
    text = path.read_text()
    assert path.name == "report.html" and text.startswith("<!doctype html>")
    # nothing is fetched from anywhere: no scripts, no external images, fonts or styles
    assert "<script" not in text and not re.search(r'(src|href)="https?://', text)
    assert full.run_id in text and full.result_checksum in text
    # user text is escaped
    assert (
        "A &lt;b&gt;scenario&lt;/b&gt; &amp; its report." in text and "<b>scenario</b>" not in text
    )
    # validity comes before results, and the warnings are spelled out
    assert text.index("Outside the validated range") < text.index("<h2>Results</h2>")
    for warning in full.metadata["validity"]["warnings"]:
        assert warning.split(":")[1].strip()[:40] in text
    for section in (
        "Uncertainty and numerical error",
        "Models, assumptions and limits",
        "Validation status",
        "Provenance",
    ):
        assert f"<h2>{section}</h2>" in text
    assert "Model-form alternative" in text and "constant albedo" in text
    assert "TH-V2" in text and "ILL-V2" in text and "pgda78.site04.surf" in text
    # seven charts, drawn with theme variables and no placeholder colour left behind
    assert text.count("<svg") == 7 and text.count('role="img"') == 7
    assert "var(--series-1)" in text and "var(--series-2)" in text and "var(--ramp-4)" in text
    assert not re.search(r"#fe[0-9a-f]{2}fe", text)
    for role, (light, dark) in TOKENS.items():
        assert f"--{role}:{light};" in text and f"--{role}:{dark};" in text
    assert "prefers-color-scheme:dark" in text


def test_report_without_thermal_or_uncertainty(runs):
    _, bare = runs
    text = write_report(bare.directory).read_text()
    assert "thermal not run" in text and "Surface temperature" not in text
    assert "No thermal ensemble was computed" in text
    assert text.count("<svg") == 4  # horizon, illumination, Earth, map
    assert "No warnings" in text and "Outside the validated range" not in text
    rows = {row[0]: row for row in key_values(bare.metadata)}
    assert rows["Visible solar disc, mean"][2:5] == ("not computed",) * 3
    assert not any("temperature" in name.lower() for name in rows)


def test_markdown_report_and_command_line(runs, capsys, tmp_path):
    full, _ = runs
    assert main(["report", str(full.directory), "--format", "md"]) == 0
    path = full.directory / "report.md"
    assert capsys.readouterr().out.strip() == str(path)
    text = path.read_text()
    assert (
        text.startswith("# report-test") and "| Quantity | Nominal |" in text and "WARNING" in text
    )
    images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    assert len(images) == 7 and all(
        (full.directory / image).stat().st_size > 5000 for image in images
    )
    elsewhere = tmp_path / "copy.html"
    assert main(["report", str(full.directory), "-o", str(elsewhere)]) == 0 and elsewhere.exists()
    assert main(["report", str(tmp_path / "not-a-run")]) == 2
    assert "not a run directory" in capsys.readouterr().err
