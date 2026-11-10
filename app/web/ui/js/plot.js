// Plotly loading, the series colours and line styles, and the shared look of
// every plot. Colours come from theme.css.

const PLOTLY_CDN = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js";

function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = src;
    s.onload = resolve;
    s.onerror = reject;
    document.head.appendChild(s);
  });
}

// ui/vendor/plotly.min.js first (offline; `serve.py --fetch-plotly`), then the CDN.
export async function loadPlotly() {
  if (window.Plotly) return window.Plotly;
  try {
    await loadScript("vendor/plotly.min.js");
  } catch {
    await loadScript(PLOTLY_CDN);
  }
  return window.Plotly;
}

export const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// fn's result with the light theme's colours, for figures on white paper.
// Synchronous, so the page never paints in between.
export function inLightTheme(fn) {
  const root = document.documentElement, theme = root.dataset.theme;
  root.dataset.theme = "light";
  try { return fn(); } finally { root.dataset.theme = theme; }
}

export const seriesColor = (i) => css(`--s${(i % 8) + 1}`);
export const DASHES = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"];
const DASH_SVG = ["", "5 3", "1.5 3", "5 3 1.5 3", "9 3", "9 3 1.5 3"];

// Small SVG of a line style, for panel rows and legends.
export function strokeSvg(i, color = "currentColor", width = 22) {
  return `<svg class="stroke" width="${width}" height="8" aria-hidden="true"><line x1="1" y1="4" x2="${width - 1}" y2="4" stroke="${color}" stroke-width="2" stroke-dasharray="${DASH_SVG[i % DASH_SVG.length]}"/></svg>`;
}

export function axisStyle() {
  return {
    gridcolor: css("--grid"),
    linecolor: css("--line-strong"),
    zeroline: false,
    showline: true,
    mirror: true,
    ticks: "outside",
    ticklen: 4,
    tickcolor: css("--line-strong"),
    tickfont: { size: 11, color: css("--muted") },
    title: { font: { size: 12.5, color: css("--ink-2") } },
  };
}

export function baseLayout() {
  return {
    font: { family: css("--font-sans"), size: 12, color: css("--ink") },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: css("--bg-card"),
    margin: { l: 62, r: 14, t: 26, b: 48 },
    hovermode: "closest",
    hoverlabel: { font: { family: css("--font-mono"), size: 11.5 }, bgcolor: css("--bg-popover"), bordercolor: css("--line-strong") },
    legend: {
      bgcolor: css("--bg-card"), bordercolor: css("--line"), borderwidth: 1,
      x: 0.98, xanchor: "right", y: 0.98, yanchor: "top", font: { size: 11.5 },
    },
    modebar: { bgcolor: "rgba(0,0,0,0)", color: css("--muted"), activecolor: css("--accent") },
    showlegend: false,
  };
}

export function logAxis(title) {
  return {
    ...axisStyle(), type: "log", exponentformat: "power", dtick: 1,
    title: { ...axisStyle().title, text: title, standoff: 6 },
  };
}

export function linAxis(title) {
  return { ...axisStyle(), title: { ...axisStyle().title, text: title, standoff: 6 } };
}

export const CONFIG = {
  displaylogo: false,
  responsive: true,
  displayModeBar: "hover",
  modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines", "zoomIn2d", "zoomOut2d"],
  toImageButtonOptions: { format: "png", scale: 2, filename: "sigmatau" },
};

// Low to high σ: low values recede into the card, high values stand out.
export const heatScale = () => [[0, css("--heat-lo")], [0.55, css("--heat-mid")], [1, css("--heat-hi")]];

export function hexAlpha(hex, a) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

// Pick an SI prefix for seconds so tick labels read like "20 ns".
export function timeUnit(maxAbs) {
  const units = [[1, "s"], [1e-3, "ms"], [1e-6, "µs"], [1e-9, "ns"], [1e-12, "ps"], [1e-15, "fs"]];
  for (const [f, u] of units) if (maxAbs >= f) return [f, u];
  return [1e-15, "fs"];
}
