"""Automatic report of a run: one self-contained HTML file, or Markdown with images.

The report states what was computed, whether it is inside the validated range, the
nominal values next to their uncertainty, the model-form alternative, the assumptions,
the validation status and the provenance. It adds no numbers of its own: everything is
read from ``results.nc`` and ``metadata.json`` of the run directory.

Design follows ``DESIGN.md``. Charts are drawn with matplotlib and embedded as inline
SVG whose colours are CSS variables, so that they follow the reader's light or dark
theme; in Markdown they are PNG files in the light theme.
"""

from __future__ import annotations

import html
import io
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from lunarenv.runner import load_run

# role -> (light, dark). Chart colours validated with the dataviz palette checker.
TOKENS: dict[str, tuple[str, str]] = {
    "surface": ("#fcfcfb", "#1a1a19"),
    "surface-2": ("#f3f3f1", "#232322"),
    "border": ("#deddd8", "#3a3a37"),
    "text": ("#0b0b0b", "#f0efec"),
    "text-2": ("#52514e", "#c3c2b7"),
    "series-1": ("#2a78d6", "#3987e5"),
    "series-2": ("#eb6834", "#d95926"),
    "ramp-1": ("#86b6ef", "#9ec5f4"),
    "ramp-2": ("#3987e5", "#5598e7"),
    "ramp-3": ("#1c5cab", "#2a78d6"),
    "ramp-4": ("#0d366b", "#184f95"),
    "neutral": ("#c9c8c2", "#4a4a46"),
    "good": ("#0ca30c", "#0ca30c"),
    "warning": ("#fab219", "#fab219"),
    "serious": ("#ec835a", "#ec835a"),
    "critical": ("#d03b3b", "#d03b3b"),
}
# Unique placeholder colours drawn by matplotlib and swapped for CSS variables in the SVG.
# The pattern #feNNfe cannot collide with a token value (a collision would repaint real ink).
_SENTINEL = {role: f"#fe{i:02x}fe" for i, role in enumerate(TOKENS, start=1)}
assert not set(_SENTINEL.values()) & {v for pair in TOKENS.values() for v in pair}
_STATE = {
    "PASS": ("good", "✓"),
    "FAIL": ("critical", "✕"),
    "PARTIAL": ("warning", "◐"),
    "PENDING": ("text-2", "…"),
    "INFO": ("text-2", "i"),
    "NOT RUN": ("text-2", "-"),
}


_DASH = "-"  # placeholder for "not applicable" in tables


def _color(role: str, themed: bool) -> str:
    return _SENTINEL[role] if themed else TOKENS[role][0]


# ------------------------------------------------------------------------------ charts
def _axes(themed: bool, height: float = 3.0) -> tuple[Figure, Any]:
    fig, ax = plt.subplots(figsize=(9.6, height))
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(_color("border", themed))
    ax.tick_params(colors=_color("text-2", themed), labelsize=9, length=3)
    ax.grid(axis="y", color=_color("border", themed), linewidth=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def _label(ax: Any, ylabel: str, themed: bool, xlabel: str = "") -> None:
    ax.set_ylabel(ylabel, color=_color("text-2", themed), fontsize=9)
    if xlabel:
        ax.set_xlabel(xlabel, color=_color("text-2", themed), fontsize=9)


def _time_axis(ax: Any) -> None:
    locator = mdates.AutoDateLocator(minticks=4, maxticks=12)  # type: ignore[no-untyped-call]
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))  # type: ignore[no-untyped-call]
    ax.margins(x=0.0)


def _legend(ax: Any, themed: bool, **kwargs: Any) -> None:
    legend = ax.legend(frameon=False, fontsize=9, **kwargs)
    for text in legend.get_texts():
        text.set_color(_color("text", themed))


def _daily_mean(time: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    days = time.astype("datetime64[D]")
    unique, inverse = np.unique(days, return_inverse=True)
    sums = np.bincount(inverse, weights=values)
    return unique, sums / np.bincount(inverse)


def build_figures(
    dataset: xr.Dataset, metadata: dict[str, Any], *, themed: bool
) -> dict[str, Figure]:
    """Charts of a run, keyed by a stable name. ``themed`` draws with placeholder colours."""

    def c(role: str) -> str:
        return _color(role, themed)

    time = dataset["time"].values
    figures: dict[str, Figure] = {}

    # horizon and Sun positions
    fig, ax = _axes(themed, 3.2)
    azimuth = dataset["azimuth"].values
    horizon = dataset["horizon_elevation"].values
    floor = min(-2.5, float(horizon.min()) - 0.3)
    ax.fill_between(azimuth, floor, horizon, color=c("neutral"), linewidth=0, label="terrain")
    visible = dataset["sun_visible_fraction"].values > 0.0
    sun_az, sun_el = dataset["sun_azimuth"].values, dataset["sun_elevation"].values
    # one polyline per visibility state, broken where the state changes or the azimuth wraps
    wrapped = np.concatenate(([False], np.abs(np.diff(sun_az)) > 180.0))
    for state, role, label in (
        (~visible, "text-2", "Sun hidden"),
        (visible, "series-1", "Sun at least partly visible"),
    ):
        elevation = np.where(state & ~wrapped, sun_el, np.nan)
        ax.plot(sun_az, elevation, color=c(role), linewidth=1.2, label=label)
    ax.set_xlim(0, 360)
    ax.set_ylim(floor, max(3.0, float(sun_el.max()) + 0.5))
    ax.set_xticks(np.arange(0, 361, 45))
    _label(ax, "elevation (deg)", themed, "grid azimuth (deg, from map +y towards +x)")
    _legend(ax, themed, loc="upper right", ncols=3)
    figures["horizon"] = fig

    # illumination through the period
    fig, ax = _axes(themed, 2.6)
    days, mean = _daily_mean(time, dataset["sun_visible_fraction"].values)
    ax.fill_between(days, 0, 100 * mean, color=c("series-1"), alpha=0.18, linewidth=0)
    ax.plot(days, 100 * mean, color=c("series-1"), linewidth=1.5)
    ax.set_ylim(0, 100)
    _time_axis(ax)
    _label(ax, "visible solar disc, daily mean (%)", themed)
    figures["illumination"] = fig

    if "earth_visible_fraction" in dataset:
        fig, ax = _axes(themed, 2.4)
        days, mean = _daily_mean(time, dataset["earth_visible_fraction"].values)
        ax.fill_between(days, 0, 100 * mean, color=c("series-1"), alpha=0.18, linewidth=0)
        ax.plot(days, 100 * mean, color=c("series-1"), linewidth=1.5)
        ax.set_ylim(0, 100)
        _time_axis(ax)
        _label(ax, "visible Earth disc, daily mean (%)", themed)
        figures["earth"] = fig

    if "surface_temperature" in dataset:
        fig, ax = _axes(themed, 2.4)
        ax.plot(time, dataset["absorbed_solar_flux"].values, color=c("series-1"), linewidth=1.2)
        ax.set_ylim(bottom=0)
        _time_axis(ax)
        _label(ax, "absorbed solar flux (W m⁻²)", themed)
        figures["flux"] = fig

        fig, ax = _axes(themed, 3.4)
        if "surface_temperature_envelope_quantiles" in dataset:
            envelope = dataset["surface_temperature_envelope_quantiles"].values  # case, q, time
            ax.fill_between(
                time,
                envelope[:, 0, :].min(axis=0),
                envelope[:, -1, :].max(axis=0),
                color=c("neutral"),
                alpha=0.6,
                linewidth=0,
                label="bounds: albedo envelope x ensemble 5-95 %",
            )
        if "surface_temperature_quantiles" in dataset:
            q = dataset["surface_temperature_quantiles"].values
            ax.fill_between(
                time,
                q[0],
                q[-1],
                color=c("series-1"),
                alpha=0.2,
                linewidth=0,
                label="5-95 % of the uncertainty ensemble",
            )
        ax.plot(
            time,
            dataset["surface_temperature"].values,
            color=c("series-1"),
            linewidth=1.5,
            label="nominal (published albedo law)",
        )
        if "surface_temperature_constant_albedo" in dataset:
            ax.plot(
                time,
                dataset["surface_temperature_constant_albedo"].values,
                color=c("series-2"),
                linewidth=1.5,
                label="alternative: constant albedo",
            )
        _time_axis(ax)
        _label(ax, "surface temperature (K)", themed)
        _legend(ax, themed, loc="upper center", ncols=3, bbox_to_anchor=(0.5, 1.16))
        figures["surface_temperature"] = fig

        fig, ax = _axes(themed, 3.0)
        depth = dataset["depth"].values
        targets = [d for d in (0.05, 0.1, 0.3, 1.0) if d <= depth[-1]]
        for i, target in enumerate(targets):
            index = int(np.argmin(np.abs(depth - target)))
            ax.plot(
                time,
                dataset["temperature"].values[:, index],
                color=c(f"ramp-{i + 1}"),
                linewidth=1.5,
                label=f"{depth[index]:.2f} m",
            )
        _time_axis(ax)
        _label(ax, "temperature at depth (K)", themed)
        _legend(ax, themed, loc="upper center", ncols=4, bbox_to_anchor=(0.5, 1.16), title=None)
        figures["depth_temperature"] = fig

    if "map_mean_visible_fraction" in dataset:
        fig, ax = _axes(themed, 4.6)
        ramp = LinearSegmentedColormap.from_list(
            "ramp", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
        )
        x, y = dataset["x"].values / 1e3, dataset["y"].values / 1e3
        image = ax.imshow(
            100 * dataset["map_mean_visible_fraction"].values,
            extent=(x[0], x[-1], y[-1], y[0]),
            cmap=ramp,
            vmin=0,
            vmax=100,
            origin="upper",
        )
        ax.grid(False)
        ax.plot(
            metadata["site"]["x_m"] / 1e3,
            metadata["site"]["y_m"] / 1e3,
            marker="+",
            color=c("text"),
            markersize=12,
        )
        bar = fig.colorbar(image, ax=ax, shrink=0.8)
        bar.ax.tick_params(colors=c("text-2"), labelsize=9)
        bar.set_label("mean visible solar disc (%)", color=c("text-2"), fontsize=9)
        bar.outline.set_edgecolor(c("border"))
        _label(ax, "polar stereographic y (km)", themed, "polar stereographic x (km)")
        figures["map"] = fig
    for figure in figures.values():
        figure.tight_layout()
    return figures


def _svg(figure: Figure, label: str) -> str:
    buffer = io.StringIO()
    with matplotlib.rc_context({"svg.fonttype": "none", "svg.hashsalt": "lunarenv"}):
        figure.savefig(buffer, format="svg", metadata={"Date": None, "Creator": None})
    plt.close(figure)
    text = buffer.getvalue()
    text = text[text.index("<svg") :]
    for role, sentinel in _SENTINEL.items():
        text = text.replace(sentinel, f"var(--{role})")
    # matplotlib writes black text; let it inherit the theme
    text = text.replace("fill: #000000", "fill: var(--text)").replace(
        'fill="#000000"', 'fill="var(--text)"'
    )
    return text.replace("<svg ", f'<svg role="img" aria-label="{html.escape(label)}" ', 1)


# -------------------------------------------------------------------------------- text
def _fmt(value: Any, digits: int = 1) -> str:
    if value is None:
        return _DASH
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def key_values(metadata: dict[str, Any]) -> list[tuple[str, str, str, str, str, str]]:
    """Rows (quantity, nominal, 5 %, median, 95 %, unit) of the headline table."""
    rows = []
    summary, uncertainty = metadata["summary"], metadata["uncertainty"]
    ill = summary["illumination"]
    ensemble = uncertainty.get("illumination")

    def spread(
        block: dict[str, Any] | None, key: str, scale: float = 1.0, digits: int = 1
    ) -> tuple[str, str, str]:
        if not block or key not in block:
            return ("not computed",) * 3
        q = block[key]["quantiles"]
        return (_fmt(scale * q[0], digits), _fmt(scale * q[1], digits), _fmt(scale * q[2], digits))

    rows.append(
        (
            "Visible solar disc, mean",
            _fmt(100 * ill["mean_visible_fraction"], 2),
            *spread(ensemble, "mean_visible_fraction", 100.0, 2),
            "%",
        )
    )
    rows.append(
        (
            "Time with any part of the Sun visible",
            _fmt(100 * ill["fraction_of_time_lit"], 1),
            _DASH,
            _DASH,
            _DASH,
            "%",
        )
    )
    rows.append(
        (
            "Longest continuous shadow",
            _fmt(ill["longest_shadow_hours"], 0),
            *spread(ensemble, "longest_shadow_hours", 1.0, 0),
            "h",
        )
    )
    rows.append(
        (
            "Longest continuous illumination",
            _fmt(ill["longest_illumination_hours"], 0),
            _DASH,
            _DASH,
            _DASH,
            "h",
        )
    )
    if "earth_fraction_of_time_visible" in ill:
        rows.append(
            (
                "Time with any part of the Earth visible",
                _fmt(100 * ill["earth_fraction_of_time_visible"], 1),
                _DASH,
                _DASH,
                _DASH,
                "%",
            )
        )
    thermal = summary.get("thermal")
    if thermal:
        spread_t = uncertainty.get("thermal")
        rows.append(
            (
                "Surface temperature, minimum",
                _fmt(thermal["surface_min_k"]),
                *spread(spread_t, "surface_min_k"),
                "K",
            )
        )
        rows.append(
            (
                "Surface temperature, maximum",
                _fmt(thermal["surface_max_k"]),
                *spread(spread_t, "surface_max_k"),
                "K",
            )
        )
        rows.append(
            (
                "Surface temperature, mean",
                _fmt(thermal["surface_mean_k"]),
                *spread(spread_t, "surface_mean_k"),
                "K",
            )
        )
        deep = thermal["at_depth"].get("1 m")
        if deep:
            rows.append(
                (
                    "Temperature at 1 m, mean",
                    _fmt(deep["mean_k"]),
                    *spread(spread_t, "mean_at_1m_k"),
                    "K",
                )
            )
        rows.append(
            (
                "Mean amplitude of a lunation cycle",
                _fmt(thermal["lunation_amplitude_mean_k"]),
                _DASH,
                _DASH,
                _DASH,
                "K",
            )
        )
        rows.append(
            (
                "Fastest heating / cooling",
                f"{_fmt(thermal['max_heating_rate_k_per_h'])} / {_fmt(thermal['max_cooling_rate_k_per_h'])}",
                _DASH,
                _DASH,
                _DASH,
                "K/h",
            )
        )
    return rows


CAPTIONS = {
    "horizon": (
        "Skyline and Sun",
        "Terrain horizon seen from the observer, with the Sun's centre at every time step. The Sun is hidden where its marks fall inside the terrain silhouette.",
    ),
    "illumination": (
        "Solar illumination",
        "Daily mean of the visible fraction of the solar disc. Gaps to 100 % are terrain shadow or polar night.",
    ),
    "earth": (
        "Earth visibility",
        "Daily mean of the visible fraction of the Earth's disc, relevant for direct radio links.",
    ),
    "flux": (
        "Absorbed sunlight",
        "Direct solar flux absorbed by the surface element on its own slope. It is zero when the Sun is below the local surface plane, even if visible over a depressed horizon.",
    ),
    "surface_temperature": (
        "Surface temperature",
        "Nominal result (Hayne 2017 albedo law) with its ensemble band, the wider bounds over the published albedo laws (cold: Bürger 2024, hot: Foote 2020) and the constant-albedo alternative. The bounds, not the nominal line, are the result to design with.",
    ),
    "depth_temperature": (
        "Temperature at depth",
        "The seasonal wave is damped and delayed with depth; at 1 m it is nearly flat.",
    ),
    "map": (
        "Illumination map",
        "Mean visible fraction of the solar disc around the site at the observer height; the cross marks the site.",
    ),
}

_CSS = """
:root{color-scheme:light dark;%(light)s}
@media (prefers-color-scheme:dark){:root{%(dark)s}}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--text);font:16px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1040px;margin:0 auto;padding:32px 16px 48px}
h1{font-size:28px;line-height:1.2;margin:0 0 8px}
h2{font-size:20px;margin:48px 0 12px;padding-top:16px;border-top:1px solid var(--border)}
h3{font-size:16px;margin:24px 0 4px}
p{margin:0 0 12px;max-width:75ch}
.caption,.meta{color:var(--text-2);font-size:13px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px}
a{color:var(--series-1)}
dl.identity{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px;margin:16px 0;font-size:14px}
dl.identity dt{color:var(--text-2)}
dl.identity dd{margin:0;overflow-wrap:anywhere}
.banner{background:var(--surface-2);border-left:4px solid var(--state);border-radius:4px;padding:12px 16px;margin:24px 0}
.banner ul{margin:8px 0 0;padding-left:20px;font-size:14px}
.badge{display:inline-flex;align-items:center;gap:4px;font-size:13px;font-weight:600;white-space:nowrap}
.badge i{font-style:normal;display:inline-grid;place-items:center;width:18px;height:18px;border-radius:4px;background:var(--state);color:var(--surface);font-size:12px}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%%;font-size:14px;margin:8px 0 16px}
th,td{text-align:left;padding:6px 12px 6px 0;border-bottom:1px solid var(--border);vertical-align:top}
th{color:var(--text-2);font-weight:600;font-size:13px;background:var(--surface)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
figure{margin:24px 0}
figure svg{width:100%%;height:auto;display:block}
details{margin:8px 0}
summary{cursor:pointer;font-weight:600;font-size:14px;padding:4px 0}
summary:focus-visible,a:focus-visible{outline:2px solid var(--series-1);outline-offset:2px}
pre{background:var(--surface-2);padding:12px;border-radius:4px;overflow-x:auto;font-size:12px}
ul.tight{margin:4px 0 12px;padding-left:20px;font-size:14px}
@media print{main{max-width:none}h2{break-after:avoid}figure{break-inside:avoid}}
"""


def _badge(status: str) -> str:
    role, glyph = _STATE.get(status, ("text-2", "?"))
    return f'<span class="badge" style="--state:var(--{role})"><i aria-hidden="true">{glyph}</i>{html.escape(status)}</span>'


def _maybe(value: float | None) -> str:
    return "not computed" if value is None else f"{value:.2g}"


def _table(headers: list[tuple[str, bool]], rows: list[list[str]]) -> str:
    head = "".join(
        f'<th scope="col"{" class=num" if numeric else ""}>{html.escape(h)}</th>'
        for h, numeric in headers
    )
    body = "".join(
        "<tr>"
        + "".join(
            f"<td{' class=num' if headers[i][1] else ''}>{cell}</td>" for i, cell in enumerate(row)
        )
        + "</tr>"
        for row in rows
    )
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def render_html(dataset: xr.Dataset, metadata: dict[str, Any]) -> str:
    """Self-contained HTML report of a run."""
    e = html.escape
    scenario, site = metadata["scenario"], metadata["site"]
    warnings = metadata["validity"]["warnings"]
    thermal_on = "thermal" in metadata["summary"]
    figures = {
        name: _svg(figure, CAPTIONS[name][0])
        for name, figure in build_figures(dataset, metadata, themed=True).items()
    }
    css = _CSS % {
        "light": "".join(f"--{role}:{light};" for role, (light, _) in TOKENS.items()),
        "dark": "".join(f"--{role}:{dark};" for role, (_, dark) in TOKENS.items()),
    }
    parts = [
        f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>lunarenv run {e(metadata['run_id'])}: {e(scenario['name'])}</title><style>{css}</style></head><body><main>",
        f"<h1>{e(scenario['name'])}</h1>",
        f'<p class="meta">lunarenv {e(metadata["lunarenv_version"])} run report, generated from run <span class="mono">{e(metadata["run_id"])}</span></p>',
    ]
    if scenario.get("description"):
        parts.append(f"<p>{e(scenario['description'])}</p>")
    period = scenario["period"]
    identity = [
        (
            "Site",
            f"{site['latitude_deg']:.4f}°, {site['longitude_deg']:.4f}° E; map ({site['x_m']:.1f}, {site['y_m']:.1f}) m; terrain {site['terrain_height_m']:.1f} m above the {site['reference_radius_m'] / 1e3:.1f} km sphere",
        ),
        ("Terrain", e(site["terrain_description"])),
        ("Observer", f"{scenario['site']['height_above_surface_m']:g} m above the terrain"),
        (
            "Period",
            f"{period['start_utc']} to {period['end_utc']}, every {period['step_minutes']} min",
        ),
        (
            "Models",
            ", ".join(f"{m['name']} {m['version']}" for m in metadata["models"])
            + ("" if thermal_on else " (thermal not run)"),
        ),
        ("Run id", f'<span class="mono">{e(metadata["run_id"])}</span>'),
        ("Result checksum", f'<span class="mono">{e(metadata["result_checksum_sha256"])}</span>'),
        ("Created", e(metadata["created_utc"])),
    ]
    parts.append(
        '<dl class="identity">'
        + "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in identity)
        + "</dl>"
    )

    if warnings:
        parts.append(
            '<div class="banner" style="--state:var(--warning)" role="note">'
            + _badge("PARTIAL").replace("PARTIAL", "Outside the validated range")
            + "<ul>"
            + "".join(f"<li>{e(w)}</li>" for w in warnings)
            + "</ul></div>"
        )
    else:
        parts.append(
            '<div class="banner" style="--state:var(--good)" role="note">'
            + _badge("PASS").replace("PASS", "No warnings")
            + "<p class=caption>Every input is inside the range the models have been validated for.</p></div>"
        )

    parts.append("<h2>Results</h2>")
    parts.append(
        "<p>Nominal values with the 5 %, 50 % and 95 % quantiles of the uncertainty ensemble where one was computed.</p>"
    )
    parts.append(
        _table(
            [
                ("Quantity", False),
                ("Nominal", True),
                ("5 %", True),
                ("Median", True),
                ("95 %", True),
                ("Unit", False),
            ],
            [[e(r[0]), *map(e, r[1:])] for r in key_values(metadata)],
        )
    )
    alternative = metadata["model_form_sensitivity"].get("constant_albedo")
    if alternative:
        parts.append(
            f"<p><strong>Model-form alternative.</strong> With a constant albedo the surface mean is {alternative['surface_mean_k']:.1f} K "
            f"({alternative['difference_surface_mean_k']:+.1f} K) and the maximum {alternative['surface_max_k']:.1f} K "
            f"({alternative['difference_surface_max_k']:+.1f} K). {e(alternative['status'].capitalize())}. This difference is not part of the ensemble quantiles above.</p>"
        )
    bounds = metadata["uncertainty"].get("thermal_envelope", {}).get("bounds")
    if bounds:
        labels = {
            "surface_min_k": "Surface minimum",
            "surface_max_k": "Surface maximum",
            "surface_mean_k": "Surface mean",
            "mean_at_1m_k": "Mean at 1 m",
        }
        parts.append("<h3>Thermal bounds</h3>")
        parts.append(
            "<p>Lower end: 5 % quantile of the coldest albedo law; upper end: 95 % quantile of the "
            "hottest (Bürger 2024 and Foote 2020, highland laws). Validation case TH-V3c tests "
            "whether Diviner polar temperatures at the site fall inside. Terrain radiation (single "
            "bounce) and Earthshine are included when the run says so; multiple scattering is "
            "not.</p>"
        )
        parts.append(
            _table(
                [("Quantity", False), ("Lower", True), ("Upper", True), ("Unit", False)],
                [
                    [e(labels[k]), f"{v[0]:.1f}", f"{v[1]:.1f}", "K"]
                    for k, v in bounds.items()
                    if k in labels
                ],
            )
        )
    grid = metadata["summary"].get("thermal", {}).get("slope_grid")
    if grid:
        lo, hi = grid["range_surface_max_k"]
        parts.append(
            f"<p><strong>Slope grid.</strong> Over the ground slopes and aspects requested and the "
            f"albedo laws, the surface maximum ranges from {lo:.1f} to {hi:.1f} K (nominal terrain). "
            f"{e(grid['note'])}</p>"
        )
    diffuse = metadata["summary"].get("thermal", {}).get("diffuse")
    if diffuse:
        case = diffuse["cases"]["nominal"]
        parts.append(
            "<p><strong>Terrain radiation and Earthshine.</strong> At the surface minimum "
            f"({case['surface_min_k']:.1f} K, {e(case['time_of_minimum_utc'])}) the terrain sends "
            f"{case['terrain_irradiance_at_minimum_w_m2']:.2g} W m⁻² and Earthshine supplies "
            f"{case['earthshine_absorbed_at_minimum_w_m2']:.2g} W m⁻² absorbed; they raise the "
            f"minimum by {case['terrain_effect_k']['surface_min']:+.2f} K and "
            f"{case['earthshine_effect_k']['surface_min']:+.2f} K. {e(diffuse['method'])}.</p>"
        )
    surfaces = metadata["summary"].get("oriented_surfaces")
    if surfaces:
        parts.append("<h3>Asset surfaces</h3>")
        parts.append(f"<p>{e(surfaces['note'])}</p>")
        parts.append(
            _table(
                [
                    ("Surface", False),
                    ("Tilt (deg)", True),
                    ("Facing (deg)", True),
                    ("Mean direct (W m⁻²)", True),
                    ("Max direct (W m⁻²)", True),
                    ("Sunlit (%)", True),
                    ("Longest dark (h)", True),
                    ("Sky view", True),
                    ("Ground view", True),
                    ("Mean terrain scattered (W m⁻²)", True),
                    ("Mean terrain IR (W m⁻²)", True),
                    ("Mean Earthshine (W m⁻²)", True),
                ],
                [
                    [
                        e(name),
                        f"{s['tilt_deg']:g}",
                        f"{s['facing_grid_azimuth_deg']:g}",
                        f"{s['mean_direct_irradiance_w_m2']:.0f}",
                        f"{s['max_direct_irradiance_w_m2']:.0f}",
                        f"{100 * s['fraction_of_time_with_direct_sun']:.0f}",
                        f"{s['longest_without_direct_sun_hours']:.0f}",
                        f"{s['sky_view_factor']:.3f}",
                        f"{s['ground_view_factor']:.3f}",
                        _maybe(s.get("mean_terrain_scattered_w_m2")),
                        _maybe(s.get("mean_terrain_infrared_w_m2")),
                        _maybe(
                            None
                            if "mean_earthshine_reflected_w_m2" not in s
                            else s["mean_earthshine_reflected_w_m2"]
                            + s["mean_earthshine_infrared_w_m2"]
                        ),
                    ]
                    for name, s in surfaces["surfaces"].items()
                ],
            )
        )
    for name, svg in figures.items():
        title, caption = CAPTIONS[name]
        parts.append(f"<figure><h3>{e(title)}</h3><p class=caption>{e(caption)}</p>{svg}</figure>")

    parts.append("<h2>Uncertainty and numerical error</h2>")
    unc = metadata["uncertainty"]
    notes = [
        f"Terrain: {unc['terrain_clones_used']} of {unc['terrain_clones_requested']} requested DEM error clones were used."
    ]
    if "illumination" in unc:
        notes.append(e(unc["illumination"]["source"]) + ".")
    if "thermal" in unc:
        notes.append(
            f"Thermal: {unc['thermal']['members']} members. {e(unc['thermal']['source'])}."
        )
    else:
        notes.append("No thermal ensemble was computed for this run.")
    parts.append('<ul class="tight">' + "".join(f"<li>{n}</li>" for n in notes) + "</ul>")
    if "thermal" in unc:
        parts.append(
            "<details><summary>Parameter intervals sampled</summary>"
            + _table(
                [("Parameter", False), ("Low", True), ("High", True), ("Basis", False)],
                [
                    [f"<code>{e(k)}</code>", f"{v['low']:g}", f"{v['high']:g}", e(v["basis"])]
                    for k, v in unc["thermal"]["parameter_intervals"].items()
                ],
            )
            + "</details>"
        )
    error_rows = [
        [e(model), f"<code>{e(k)}</code>", e(str(v))]
        for model, block in metadata["numerical_error"].items()
        for k, v in block.items()
    ]
    parts.append(
        "<h3>Numerical error estimates</h3><p class=caption>Estimated from convergence studies at the default settings; reported separately from input uncertainty.</p>"
        + _table([("Model", False), ("Quantity", False), ("Estimate", False)], error_rows)
    )

    parts.append("<h2>Models, assumptions and limits</h2>")
    for model in metadata["models"]:
        parts.append(
            f"<h3>{e(model['name'])} {e(model['version'])}</h3><p>{e(model['summary'])}</p>"
        )
        for heading in ("assumptions", "simplifications", "validity_range", "uncertainty_sources"):
            parts.append(
                f"<details{' open' if heading == 'simplifications' else ''}><summary>{heading.replace('_', ' ').capitalize()}</summary><ul class=tight>"
                + "".join(f"<li>{e(i)}</li>" for i in model[heading])
                + "</ul></details>"
            )

    parts.append("<h2>Validation status</h2>")
    validation = metadata["validation"]
    counts = ", ".join(f"{n} {s}" for s, n in sorted(validation["counts"].items()))
    parts.append(
        f"<p>Stored results of {e(str(validation['generated_utc']))}: {e(counts)}. Sources and tolerances are in VALIDATION.md.</p>"
    )
    parts.append(
        _table(
            [("Case", False), ("Status", False), ("Result", False)],
            [
                [f"<code>{e(c['id'])}</code>", _badge(c["status"]), e(c["summary"])]
                for c in validation["cases"]
            ],
        )
    )

    parts.append("<h2>Provenance</h2>")
    parts.append(
        _table(
            [("Data", False), ("Version", False), ("SHA-256", False), ("Citation", False)],
            [
                [
                    f"<code>{e(d['id'])}</code>",
                    e(d["version"]),
                    f"<code>{e(d['sha256'][:16])}</code>",
                    e(d["citation"]),
                ]
                for d in metadata["data"]
            ],
        )
    )
    environment = metadata["environment"]
    parts.append(
        "<details><summary>Software environment</summary>"
        + _table(
            [("Component", False), ("Version", False)],
            [[e(k), e(str(v))] for k, v in environment.items()],
        )
        + "</details>"
    )
    parts.append(
        f"<details><summary>Resolved scenario</summary><pre>{e(json.dumps(scenario, indent=1))}</pre></details>"
    )
    variables = [
        [f"<code>{e(k)}</code>", e(", ".join(v["dims"])), e(v["units"]), e(v["long_name"])]
        for k, v in metadata["variables"].items()
    ]
    parts.append(
        "<details><summary>Variables in results.nc</summary>"
        + _table(
            [("Name", False), ("Dimensions", False), ("Unit", False), ("Meaning", False)], variables
        )
        + "</details>"
    )
    parts.append("</main></body></html>")
    return "\n".join(parts)


def render_markdown(dataset: xr.Dataset, metadata: dict[str, Any], figure_dir: Path) -> str:
    """Markdown report; charts are written as PNG files into ``figure_dir``."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    scenario, site = metadata["scenario"], metadata["site"]
    lines = [
        f"# {scenario['name']}",
        "",
        f"lunarenv {metadata['lunarenv_version']} run `{metadata['run_id']}`, result checksum `{metadata['result_checksum_sha256']}`.",
        "",
        f"- Site: {site['latitude_deg']:.4f} deg, {site['longitude_deg']:.4f} deg E; map ({site['x_m']:.1f}, {site['y_m']:.1f}) m; {site['terrain_description']}",
        f"- Observer: {scenario['site']['height_above_surface_m']:g} m above the terrain",
        f"- Period: {scenario['period']['start_utc']} to {scenario['period']['end_utc']}, every {scenario['period']['step_minutes']} min",
        "",
        "## Validity",
        "",
    ]
    warnings = metadata["validity"]["warnings"]
    lines += [f"- WARNING: {w}" for w in warnings] or ["No warnings."]
    lines += [
        "",
        "## Results",
        "",
        "| Quantity | Nominal | 5 % | Median | 95 % | Unit |",
        "|---|---:|---:|---:|---:|---|",
    ]
    lines += [f"| {' | '.join(row)} |" for row in key_values(metadata)]
    alternative = metadata["model_form_sensitivity"].get("constant_albedo")
    if alternative:
        lines += [
            "",
            f"Model-form alternative (constant albedo): surface mean {alternative['surface_mean_k']:.1f} K ({alternative['difference_surface_mean_k']:+.1f} K), maximum {alternative['surface_max_k']:.1f} K ({alternative['difference_surface_max_k']:+.1f} K). {alternative['status'].capitalize()}.",
        ]
    for name, figure in build_figures(dataset, metadata, themed=False).items():
        path = figure_dir / f"{name}.png"
        figure.savefig(path, dpi=130, facecolor=TOKENS["surface"][0])
        plt.close(figure)
        title, caption = CAPTIONS[name]
        lines += ["", f"### {title}", "", caption, "", f"![{title}]({figure_dir.name}/{path.name})"]
    lines += ["", "## Validation status", "", "| Case | Status | Result |", "|---|---|---|"]
    lines += [
        f"| {c['id']} | {c['status']} | {c['summary'].replace('|', '/')} |"
        for c in metadata["validation"]["cases"]
    ]
    lines += [
        "",
        "## Provenance",
        "",
        "| Data | Version | SHA-256 | Citation |",
        "|---|---|---|---|",
    ]
    lines += [
        f"| `{d['id']}` | {d['version']} | `{d['sha256'][:16]}` | {d['citation']} |"
        for d in metadata["data"]
    ]
    lines += [
        "",
        "Assumptions, simplifications and the full provenance are in `metadata.json`.",
        "",
    ]
    return "\n".join(lines)


def write_report(
    run_directory: str | Path, *, fmt: str = "html", output: str | Path | None = None
) -> Path:
    """Write the report of a run next to its results and return the path."""
    directory = Path(run_directory)
    if not (directory / "metadata.json").exists():
        raise FileNotFoundError(f"{directory} is not a run directory (no metadata.json)")
    dataset, metadata = load_run(directory)
    if fmt == "html":
        path = Path(output) if output else directory / "report.html"
        path.write_text(render_html(dataset, metadata), encoding="utf-8")
    elif fmt in ("md", "markdown"):
        path = Path(output) if output else directory / "report.md"
        path.write_text(
            render_markdown(dataset, metadata, path.parent / "report_figures"), encoding="utf-8"
        )
    else:
        raise ValueError("format must be 'html' or 'md'")
    return path
