// Dynamic page: σ(t, τ) over sliding windows. One clock: a map (or a 3D
// waterfall) and the selected cell. Several clocks: one map per clock side by
// side on a shared color scale, since maps cannot be overlaid. Maps are
// computed by the core when this page needs them and recomputed when the
// estimator, window or record changes; no Run is needed.
import { state, update, visibleClocks, multi } from "../state.js";
import { dynamic as fetchDynamic, dynamicKey } from "../core.js";
import { loadPlotly, baseLayout, linAxis, logAxis, CONFIG, seriesColor, css, heatScale } from "../plot.js";
import { sci, num, thousands, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
const inflight = new Map();     // clock id → key being fetched
const FRACS = [4, 8, 16, 32];
const NAME = { dadev: "DADEV", dhdev: "DHDEV" };

export function controls(inline = false) {
  const d = state.dyn;
  const est = `<div class="seg" role="group" aria-label="Estimator">${Object.entries(NAME).map(([k, v]) =>
    `<button type="button" data-ctl="dev" data-v="${k}" class="${d.dev === k ? "on" : ""}">${v}</button>`).join("")}</div>`;
  if (inline) {
    return `${est}<label class="tool-field">Window <select data-ctl="frac">${FRACS.map((f) =>
      `<option value="${f}" ${f === d.frac ? "selected" : ""}>1/${f} of record</option>`).join("")}</select></label>
      <span class="caption grow-right">Maps can't be overlaid, so each clock gets its own. Same color scale throughout.</span>`;
  }
  const c = visibleClocks()[0];
  const n = c?.record?.n ?? 0;
  const w = Math.max(16, Math.floor(n / d.frac));
  return `<h2>Dynamic deviation</h2>
    <div class="ctl-block"><label class="lbl">Estimator</label>${est}</div>
    <div class="ctl-block"><label class="lbl">Window length</label>
      <div class="pills">${FRACS.map((f) => `<button type="button" class="pill ${f === d.frac ? "on" : ""}" data-ctl="frac" data-v="${f}">1/${f}</button>`).join("")}</div>
      ${n ? `<div class="muted small">${thousands(w)} samples per window · ${2 * d.frac - 1} windows</div>` : ""}</div>
    <div class="ctl-block"><label class="lbl">View</label>
      <div class="seg" role="group" aria-label="View"><button type="button" data-ctl="view" data-v="map" class="${d.view === "map" ? "on" : ""}">Map</button><button type="button" data-ctl="view" data-v="3d" class="${d.view === "3d" ? "on" : ""}">3D waterfall</button></div></div>
    <div class="rule"></div>
    <p class="muted small">Changes here recompute this page right away. A clock whose noise does not change shows the same colors in every row; a change during the record shows as a change down a column.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "dev" && e.type === "click") update((s) => { s.dyn.dev = t.dataset.v; s.dyn.cell = null; });
  if (k === "view" && e.type === "click") update((s) => { s.dyn.view = t.dataset.v; });
  if (k === "frac" && (e.type === "click" && t.dataset.v || e.type === "change" && t.tagName === "SELECT")) {
    const v = Number(t.dataset.v ?? t.value);
    update((s) => { s.dyn.frac = v; s.dyn.cell = null; });
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar" data-ref="toolbar"></div>
      <div class="dyn-area" data-ref="area"></div>
    </div>
    <div class="card side-card dyn-cell" data-ref="cell" hidden></div>`;
  root.addEventListener("click", (e) => {
    if (e.target.closest("[data-act=load]")) ctx.goData();
  });
  root.addEventListener("change", (e) => {
    if (e.target.dataset.ref === "slice-tau") update((s) => { s.dyn.sliceTau = Number(e.target.value); });
  });
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for any map that is missing or out of date.
function ensureMaps(clocks) {
  for (const c of clocks) {
    if (!c.record) continue;
    const key = dynamicKey(c);
    if (c.dynamic?.key === key || inflight.get(c.id) === key) continue;
    inflight.set(c.id, key);
    fetchDynamic(c).then((out) => {
      if (inflight.get(c.id) !== key) return;
      inflight.delete(c.id);
      update(() => { c.dynamic = { key, ...out.dynamic }; c.dynamicError = null; });
    }).catch((err) => {
      if (inflight.get(c.id) !== key) return;
      inflight.delete(c.id);
      update(() => { c.dynamic = { key, error: err.message }; });
    });
  }
}

const log10 = (v) => (v > 0 ? Math.log10(v) : NaN);
const tauText = (t) => (t >= 1e4 ? `${num(t / 1000, 3)}k` : num(t, 4));
const timeText = (t) => (t >= 1e4 ? `${thousands(Math.round(t / 1000))}k` : thousands(Math.round(t)));

function heat(c, d, zmin, zmax, showscale, compact) {
  return {
    type: "heatmap",
    z: d.map.map((row) => row.map(log10)),
    x: d.tau.map((_, j) => j),
    y: d.t,
    zmin, zmax,
    colorscale: heatScale(),
    xgap: 1, ygap: 1,
    showscale,
    colorbar: showscale ? {
      thickness: 14, outlinewidth: 0, len: 0.9,
      tickvals: range(Math.ceil(zmin), Math.floor(zmax)),
      ticktext: range(Math.ceil(zmin), Math.floor(zmax)).map((k) => `10<sup>${k}</sup>`),
      tickfont: { size: 10.5, color: css("--muted") },
    } : undefined,
    customdata: d.map.map((row, i) => row.map((v, j) => [d.tau[j], sci(v, 3), d.t[i]])),
    hovertemplate: `${compact ? `<b>${h(c.name)}</b><br>` : ""}t = %{customdata[2]:.6g} s<br>τ = %{customdata[0]:.4g} s<br>${NAME[d.dev]} = %{customdata[1]}<extra></extra>`,
  };
}

const sceneAxis = (text) => ({ title: { text }, color: css("--ink-2"), gridcolor: css("--grid"), zerolinecolor: css("--line-strong"), linecolor: css("--line-strong") });

function range(a, b) { const out = []; for (let k = a; k <= b; k++) out.push(k); return out; }

function mapLayout(d, compact) {
  const ticks = d.tau.map((_, j) => j).filter((j) => j % Math.max(1, Math.ceil(d.tau.length / 6)) === 0);
  return {
    ...baseLayout(),
    margin: compact ? { l: 46, r: 8, t: 6, b: 42 } : { l: 70, r: 10, t: 10, b: 48 },
    xaxis: { ...linAxis(compact ? "τ (s)" : "Averaging time τ (s)"), tickvals: ticks, ticktext: ticks.map((j) => tauText(d.tau[j])), showgrid: false, mirror: false },
    yaxis: { ...linAxis(compact ? "" : "Window center t (s)"), autorange: "reversed", tickvals: d.t, ticktext: d.t.map(timeText), nticks: 6, showgrid: false, mirror: false },
    plot_bgcolor: css("--bg-card"),
    uirevision: `${d.dev}-${d.window}`,
  };
}

function zRange(maps) {
  let lo = Infinity, hi = -Infinity;
  for (const d of maps) for (const row of d.map) for (const v of row) {
    const l = log10(v);
    if (Number.isFinite(l)) { lo = Math.min(lo, l); hi = Math.max(hi, l); }
  }
  return lo < hi ? [lo, hi] : [lo - 0.5, hi + 0.5];
}

export function refresh() {
  if (!root) return;
  const several = multi();
  const tb = root.querySelector("[data-ref=toolbar]");
  tb.hidden = !several;
  if (several && !tb.contains(document.activeElement)) tb.innerHTML = controls(true);
  const area = root.querySelector("[data-ref=area]");
  const cell = root.querySelector("[data-ref=cell]");
  const clocks = visibleClocks().filter((c) => c.record);

  if (!clocks.length) {
    cell.hidden = true;
    area.innerHTML = `<div class="card plot-card"><div class="empty"><h2>No record yet</h2><p>Load or synthesize a record on the Data page.</p><button class="btn" data-act="load">Go to Data</button></div></div>`;
    return;
  }
  ensureMaps(clocks);
  const ready = clocks.filter((c) => c.dynamic?.key === dynamicKey(c) && !c.dynamic.error);
  const failed = clocks.filter((c) => c.dynamic?.key === dynamicKey(c) && c.dynamic.error);
  if (!Plotly) return;

  if (!several) {
    const c = clocks[0];
    const d = ready[0]?.dynamic;
    let card = area.querySelector(".dyn-single");
    if (!card) {
      area.innerHTML = `<div class="card plot-card dyn-single"><div class="card-row"><h2 data-ref="title"></h2><span class="muted small" data-ref="sub"></span></div><div class="plot" data-ref="map"></div><div class="empty" data-ref="wait" hidden></div></div>`;
      card = area.querySelector(".dyn-single");
      bindClick(card.querySelector("[data-ref=map]"));
    }
    card.querySelector("[data-ref=title]").textContent = `${NAME[state.dyn.dev]} by window`;
    card.querySelector("[data-ref=sub]").textContent = state.dyn.view === "map" ? "Hover a cell for its value · click to pin it" : "Drag to turn the waterfall";
    const wait = card.querySelector("[data-ref=wait]");
    const plot = card.querySelector("[data-ref=map]");
    wait.hidden = !!d;
    wait.innerHTML = d ? "" : `<p>${failed.length ? h(c.dynamic.error) : "Computing the map…"}</p>`;
    plot.style.visibility = d ? "visible" : "hidden";
    if (!d) { cell.hidden = true; return; }
    const [zmin, zmax] = zRange([d]);
    renderCell(cell, c, d); // lay the card out before Plotly sizes the map
    if (state.dyn.view === "3d") {
      Plotly.react(plot, [{
        type: "surface", x: d.tau.map(Math.log10), y: d.t, z: d.map.map((r) => r.map(log10)),
        colorscale: heatScale(), cmin: zmin, cmax: zmax, showscale: false,
        hovertemplate: "τ = 10^%{x:.2f} s<br>t = %{y:.6g} s<br>log₁₀σ = %{z:.2f}<extra></extra>",
      }], {
        ...baseLayout(), margin: { l: 0, r: 0, t: 0, b: 0 },
        scene: {
          xaxis: sceneAxis("log₁₀ τ (s)"), yaxis: sceneAxis("t (s)"), zaxis: sceneAxis("log₁₀ σ"),
          camera: { eye: { x: 1.6, y: -1.6, z: 0.9 } },
        },
        uirevision: "3d",
      }, CONFIG);
    } else {
      const layout = mapLayout(d, false);
      const pin = state.dyn.cell;
      layout.shapes = pin && pin.i < d.t.length && pin.j < d.tau.length ? [{
        type: "rect", xref: "x", yref: "y", x0: pin.j - 0.5, x1: pin.j + 0.5,
        y0: d.t[pin.i] - (d.step * c.record.tau0) / 2, y1: d.t[pin.i] + (d.step * c.record.tau0) / 2,
        line: { color: css("--s2"), width: 2 },
      }] : [];
      Plotly.react(plot, [heat(c, d, zmin, zmax, true, false)], layout, CONFIG);
    }
    return;
  }

  // Several clocks: one map each, shared color scale.
  cell.hidden = true;
  const [zmin, zmax] = zRange(ready.map((c) => c.dynamic));
  const ids = clocks.map((c) => c.id).join(",");
  if (area.dataset.ids !== ids) {
    area.dataset.ids = ids;
    area.innerHTML = `<div class="dyn-multi"><div class="dyn-row"><div class="dyn-grid" style="--n:${clocks.length}">${clocks.map((c) => `
      <div class="card dyn-card" data-clock="${c.id}"><h3><span class="clock-dot" style="background:${seriesColor(c.slot)}"></span>${h(c.name)}</h3>
        <div class="plot"></div><div class="empty" hidden></div></div>`).join("")}</div>
      <div class="dyn-scale"><div class="plot" data-ref="scale"></div></div></div>
      <div class="card dyn-slice"><div class="card-row"><h2>Compare at τ</h2><select data-ref="slice-tau" aria-label="τ"></select>
        <span class="muted small grow-right" data-ref="slice-sub"></span></div><div class="plot" data-ref="slice"></div></div></div>`;
    area.querySelectorAll(".dyn-card .plot").forEach(bindSliceClick);
    // Plotly only follows window resizes; follow the row too.
    const row = area.querySelector(".dyn-row");
    new ResizeObserver(() => row.querySelectorAll(".plot").forEach((p) => { if (p.layout && p.offsetParent) Plotly.Plots.resize(p); })).observe(row);
  }
  for (const c of clocks) {
    const box = area.querySelector(`.dyn-card[data-clock="${c.id}"]`);
    const plot = box.querySelector(".plot"), wait = box.querySelector(".empty");
    const ok = ready.includes(c);
    wait.hidden = ok;
    wait.innerHTML = ok ? "" : `<p>${c.dynamic?.error ? h(c.dynamic.error) : "Computing…"}</p>`;
    plot.style.visibility = ok ? "visible" : "hidden";
    if (ok) Plotly.react(plot, [heat(c, c.dynamic, zmin, zmax, false, true)], mapLayout(c.dynamic, true), { ...CONFIG, displayModeBar: false });
  }
  // A color bar on its own, shared by every map.
  if (ready.length) {
    Plotly.react(area.querySelector("[data-ref=scale]"), [{
      type: "heatmap", z: [[zmin, zmax]], zmin, zmax, colorscale: heatScale(), showscale: true, hoverinfo: "skip", opacity: 0,
      colorbar: {
        thickness: 14, outlinewidth: 0, len: 1, x: 0, xanchor: "left",
        tickvals: range(Math.ceil(zmin), Math.floor(zmax)),
        ticktext: range(Math.ceil(zmin), Math.floor(zmax)).map((k) => `10<sup>${k}</sup>`),
        tickfont: { size: 10.5, color: css("--muted") },
      },
    }], {
      ...baseLayout(), margin: { l: 0, r: 0, t: 20, b: 40 }, plot_bgcolor: "rgba(0,0,0,0)",
      xaxis: { visible: false }, yaxis: { visible: false },
    }, { ...CONFIG, displayModeBar: false, staticPlot: true });
  }
  drawSlice(area, ready);
}

// ── Several clocks: σ over time at one τ ──
// One line per clock through its windows, so a clock that drifts or steps
// stands out against the others at the same τ. The dotted lines are each
// clock's whole-record value from the last run.

const sameTau = (a, b) => Math.abs(a - b) <= 1e-9 * Math.max(Math.abs(a), Math.abs(b));

function nearestTau(taus, tau) {
  let best = -1, err = Math.log(Math.SQRT2) + 1e-9;
  taus.forEach((t, j) => { const e = Math.abs(Math.log(t / tau)); if (e <= err) { best = j; err = e; } });
  return best;
}

function sliceTaus(ready) {
  const all = [];
  for (const c of ready) for (const t of c.dynamic.tau) if (!all.some((u) => sameTau(u, t))) all.push(t);
  return all.sort((a, b) => a - b);
}

function drawSlice(area, ready) {
  const box = area.querySelector(".dyn-slice");
  if (!box) return;
  const taus = sliceTaus(ready);
  // The card keeps its place while maps compute, so the maps above it do not
  // change size when it fills in.
  box.classList.toggle("waiting", !taus.length);
  if (!taus.length) return;
  if (!taus.some((t) => state.dyn.sliceTau && sameTau(t, state.dyn.sliceTau))) state.dyn.sliceTau = taus[Math.min(2, taus.length - 1)];
  const tau = state.dyn.sliceTau;
  const sel = box.querySelector("[data-ref=slice-tau]");
  if (document.activeElement !== sel) {
    sel.innerHTML = taus.map((t) => `<option value="${t}" ${sameTau(t, tau) ? "selected" : ""}>${num(t, 4)} s</option>`).join("");
  }
  const whole = state.dyn.dev === "dhdev" ? "hdev" : "adev";
  const traces = [];
  const shapes = [];
  const missing = [];
  for (const c of ready) {
    const d = c.dynamic;
    // Clocks with another τ₀ have no exact match: take the nearest τ within
    // a half octave and name it in the legend; list the rest.
    const j = nearestTau(d.tau, tau);
    if (j < 0) { missing.push(c.name); continue; }
    const exact = sameTau(d.tau[j], tau);
    const color = seriesColor(c.slot);
    traces.push({
      x: d.t, y: d.map.map((row) => row[j]), name: exact ? c.name : `${c.name} (τ = ${num(d.tau[j], 4)} s)`, mode: "lines+markers",
      line: { color, width: 2 }, marker: { color, size: 5 },
      hovertemplate: `<b>${h(c.name)}</b><br>t = %{x:.6g} s<br>${NAME[d.dev]} = %{y:.3e}<extra></extra>`,
    });
    const r = c.result?.results.find((x) => x.dev === whole);
    const i = r ? r.tau.findIndex((t) => sameTau(t, d.tau[j])) : -1;
    if (i >= 0 && r.sigma[i] > 0) {
      shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: r.sigma[i], y1: r.sigma[i], line: { color, width: 1, dash: "dot" } });
    }
  }
  for (const c of ready) {
    const plot = area.querySelector(`.dyn-card[data-clock="${c.id}"] .plot`);
    const j = nearestTau(c.dynamic.tau, tau);
    if (plot?.layout) Plotly.relayout(plot, { shapes: j < 0 ? [] : [{ type: "rect", xref: "x", yref: "paper", x0: j - 0.5, x1: j + 0.5, y0: 0, y1: 1, line: { color: css("--accent"), width: 2 } }] });
  }
  const none = missing.length ? `No τ near ${num(tau, 4)} s for ${missing.join(", ")} · ` : "";
  box.querySelector("[data-ref=slice-sub]").textContent = `${none}${NAME[state.dyn.dev]} in each window · dotted = whole record (${whole.toUpperCase()}, last run) · click a map column to change τ`;
  Plotly.react(box.querySelector("[data-ref=slice]"), traces, {
    ...baseLayout(), margin: { l: 70, r: 14, t: 8, b: 44 },
    xaxis: linAxis("Window center t (s)"), yaxis: { ...logAxis(`${NAME[state.dyn.dev]} at τ = ${num(tau, 4)} s`), dtick: undefined },
    showlegend: true, legend: { ...baseLayout().legend, orientation: "h", x: 0, xanchor: "left", y: 1.02, yanchor: "bottom", bgcolor: "rgba(0,0,0,0)", borderwidth: 0 },
    shapes, uirevision: "slice",
  }, CONFIG);
}

let sliceBound = new WeakSet();
function bindSliceClick(plot) {
  if (sliceBound.has(plot)) return;
  sliceBound.add(plot);
  const tryBind = () => {
    if (!plot.on) return setTimeout(tryBind, 100);
    plot.on("plotly_click", (ev) => {
      const j = ev.points?.[0]?.x;
      const id = Number(plot.closest("[data-clock]")?.dataset.clock);
      const c = state.clocks.find((k) => k.id === id);
      const tau = Number.isInteger(j) ? c?.dynamic?.tau?.[j] : undefined;
      if (tau) update((s) => { s.dyn.sliceTau = tau; });
    });
  };
  tryBind();
}

let bound = new WeakSet();
function bindClick(plot) {
  if (bound.has(plot)) return;
  bound.add(plot);
  const tryBind = () => {
    if (!plot.on) return setTimeout(tryBind, 100);
    plot.on("plotly_click", (ev) => {
      const p = ev.points?.[0];
      if (!p || state.dyn.view !== "map") return;
      const i = p.pointIndex?.[0] ?? p.pointNumber?.[0], j = p.pointIndex?.[1] ?? p.pointNumber?.[1];
      if (Number.isInteger(i) && Number.isInteger(j)) update((s) => { s.dyn.cell = { i, j }; });
    });
  };
  tryBind();
}

function renderCell(box, c, d) {
  const pin = state.dyn.cell;
  // The card stays in place so the map does not resize when a cell is pinned.
  box.hidden = state.dyn.view !== "map";
  if (!pin || pin.i >= d.t.length || pin.j >= d.tau.length) {
    box.innerHTML = `<div class="card-head"><h2>Selected cell</h2></div><p class="muted small cell-kv">Click a cell to pin it and compare it with the whole record.</p>`;
    return;
  }
  const tau = d.tau[pin.j];
  const whole = c.result?.results.find((r) => r.dev === (d.dev === "dhdev" ? "hdev" : "adev"));
  const wi = whole ? whole.tau.findIndex((t) => Math.abs(t - tau) <= 1e-9 * tau) : -1;
  box.innerHTML = `<div class="card-head"><h2>Selected cell</h2></div>
    <div class="cell-kv">
      <div class="kv"><span>Window center</span><span class="mono">${thousands(Math.round(d.t[pin.i]))} s</span></div>
      <div class="kv"><span>τ</span><span class="mono">${num(tau, 5)} s</span></div>
      <div class="kv"><span>${NAME[d.dev]}</span><span class="mono">${sci(d.map[pin.i][pin.j], 3)}</span></div>
      <div class="kv" title="${d.dev === "dhdev" ? "HDEV" : "ADEV"} of the whole record at this τ, from the last run"><span>Whole record</span><span class="mono">${wi >= 0 ? sci(whole.sigma[wi], 3) : "–"}</span></div>
    </div>`;
}

export function hint() {
  if (multi()) return "Rows = window center time (s) · columns = τ · stronger color = less stable";
  return state.dyn.view === "map" ? "Click a cell to pin it · switch to 3D waterfall to compare curve shapes" : "Drag to turn · scroll to zoom";
}
