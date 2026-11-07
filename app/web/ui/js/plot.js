// Plotly loading and the shared look of every plot. Colors come from the
// --plot-* tokens in theme.css.

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

// vendor/plotly.min.js first (offline; `serve.py --fetch-plotly` puts it
// there), then the CDN.
export async function loadPlotly() {
  if (window.Plotly) return window.Plotly;
  try {
    await loadScript("vendor/plotly.min.js");
  } catch {
    await loadScript(PLOTLY_CDN);
  }
  return window.Plotly;
}

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export function palette() {
  return [1, 2, 3, 4, 5, 6, 7, 8].map((i) => css(`--plot-${i}`));
}

export const DASHES = ["solid", "dash", "dot", "dashdot", "longdash"];
export const SYMBOLS = ["circle", "square", "diamond", "triangle-up", "cross"];

export function baseLayout(extra = {}) {
  const axis = {
    gridcolor: css("--plot-grid"),
    linecolor: css("--border-strong"),
    zeroline: false,
    showline: true,
    mirror: true,
    ticks: "outside",
    tickcolor: css("--border-strong"),
    color: css("--plot-axis"),
  };
  return {
    font: { family: css("--font-sans"), size: 13, color: css("--text") },
    paper_bgcolor: css("--plot-bg"),
    plot_bgcolor: css("--plot-bg"),
    margin: { l: 72, r: 20, t: 34, b: 56 },
    hovermode: "closest",
    hoverlabel: { font: { family: css("--font-mono"), size: 12 } },
    legend: {
      bgcolor: "rgba(255,253,249,0.9)",
      bordercolor: css("--border"),
      borderwidth: 1,
      x: 0.99, xanchor: "right", y: 0.99, yanchor: "top",
    },
    xaxis: { ...axis },
    yaxis: { ...axis },
    _axis: axis,
    ...extra,
  };
}

export const CONFIG = {
  displaylogo: false,
  responsive: true,
  modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines"],
  toImageButtonOptions: { format: "png", scale: 2, filename: "sigmatau" },
};

// Log-axis tick labels as 10^n with superscripts.
export function logAxis(title) {
  return { type: "log", exponentformat: "power", dtick: 1, minor: { ticks: "outside", ticklen: 3, showgrid: false },
    title: { text: title, standoff: 8 } };
}
