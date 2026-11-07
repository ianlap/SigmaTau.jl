// Stability page: estimator picker and options on the left, σ(τ) plot in the
// middle, per-τ table on the right. Built once; refresh() updates it in place.
import { state, update, activeClocks } from "../state.js";
import { ESTIMATORS, byId, label, GRIDS, CONFIDENCES, noiseLabel } from "../estimators.js";
import { baseLayout, CONFIG, palette, DASHES, SYMBOLS, logAxis } from "../plot.js";
import { sci, num, escapeHtml, thousands } from "../format.js";

const AXIS_TITLE = {
  sigma_y: "σ<sub>y</sub>(τ)",
  sigma_x: "σ<sub>x</sub>(τ) (s)",
};

export function mountStability(root, { run }) {
  root.innerHTML = `
  <div class="stability">
    <aside class="side panel">
      <section class="side-section">
        <div class="side-head"><h2>Estimators</h2><span class="muted" data-ref="count"></span></div>
        <div class="chips" data-ref="chips"></div>
        <div class="add-wrap">
          <button class="btn btn-ghost add-btn" type="button" data-ref="add-btn" aria-expanded="false">+ Add estimator</button>
          <div class="picker" data-ref="picker" hidden></div>
        </div>
      </section>
      <section class="side-section">
        <h2>Options</h2>
        <label class="field"><span>τ grid</span>
          <select data-opt="grid">${GRIDS.map(([v, t]) => `<option value="${v}">${t}</option>`).join("")}</select></label>
        <div class="field">
          <label class="check"><input type="checkbox" data-opt="ci"> Confidence intervals</label>
          <select data-opt="confidence" aria-label="Confidence level">${CONFIDENCES.map(([v, t]) => `<option value="${v}">${t}</option>`).join("")}</select>
        </div>
        <label class="check"><input type="checkbox" data-opt="band"> Show intervals as a band</label>
        <label class="check"><input type="checkbox" data-opt="bias"> Bias-correct totals and Thêo1</label>
        <p class="hint">Detrending, scaling and τ₀ are set per clock on the <a href="#data">Data</a> page.</p>
      </section>
      <div class="side-foot">
        <button class="btn btn-primary btn-run" type="button" data-ref="run">Run analysis <kbd>F5</kbd></button>
        <p class="muted small" data-ref="last-run"></p>
      </div>
    </aside>

    <section class="center">
      <div class="toolbar">
        <div class="seg" role="radiogroup" aria-label="Plot axis" data-ref="axis-seg">
          <button type="button" data-axis="sigma_y">Frequency σ<sub>y</sub>(τ) <span class="seg-n"></span></button>
          <button type="button" data-axis="sigma_x">Time σ<sub>x</sub>(τ) <span class="seg-n"></span></button>
        </div>
        <span class="muted small" data-ref="plot-caption"></span>
        <span class="spacer"></span>
        <div class="clock-toggles" data-ref="clock-toggles"></div>
        <div class="seg" role="radiogroup" aria-label="Clock layout" data-ref="layout-seg">
          <button type="button" data-layout="overlay">Overlay</button>
          <button type="button" data-layout="side">Side by side</button>
        </div>
      </div>
      <div class="plot-card panel">
        <div class="plot" data-ref="plot"></div>
        <div class="empty" data-ref="empty" hidden></div>
      </div>
    </section>

    <aside class="results panel">
      <div class="results-head">
        <h2>Per-τ results</h2>
        <div class="results-pickers">
          <select data-ref="table-clock" aria-label="Clock"></select>
          <select data-ref="table-dev" aria-label="Estimator"></select>
        </div>
      </div>
      <div class="table-wrap"><table class="num-table" data-ref="table"></table></div>
      <p class="results-foot muted small" data-ref="table-foot"></p>
    </aside>
  </div>`;

  const $ = (ref) => root.querySelector(`[data-ref="${ref}"]`);
  const plotEl = $("plot");

  // ── Events ──
  $("chips").addEventListener("click", (e) => {
    const b = e.target.closest("[data-remove]");
    if (!b) return;
    update((s) => { s.devs = s.devs.filter((d) => d !== b.dataset.remove); }, "analysis");
  });

  const picker = $("picker"), addBtn = $("add-btn");
  const closePicker = () => { picker.hidden = true; addBtn.setAttribute("aria-expanded", "false"); };
  addBtn.addEventListener("click", () => {
    picker.hidden = !picker.hidden;
    addBtn.setAttribute("aria-expanded", String(!picker.hidden));
    if (!picker.hidden) renderPicker();
  });
  picker.addEventListener("click", (e) => {
    const b = e.target.closest("[data-add]");
    if (!b) return;
    update((s) => { s.devs = [...s.devs, b.dataset.add]; }, "analysis");
    renderPicker();
    if (!picker.querySelector("[data-add]")) closePicker();
  });
  document.addEventListener("click", (e) => {
    if (!picker.hidden && !e.target.closest(".add-wrap")) closePicker();
  });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") closePicker(); });

  root.querySelectorAll("[data-opt]").forEach((el) => {
    el.addEventListener("change", () => {
      const k = el.dataset.opt;
      const v = el.type === "checkbox" ? el.checked : k === "confidence" ? Number(el.value) : el.value;
      update((s) => { s[k] = v; }, k === "band" ? "view" : "analysis");
    });
  });

  $("run").addEventListener("click", () => run());
  $("axis-seg").addEventListener("click", (e) => {
    const b = e.target.closest("[data-axis]");
    if (b) update((s) => { s.axis = b.dataset.axis; });
  });
  $("layout-seg").addEventListener("click", (e) => {
    const b = e.target.closest("[data-layout]");
    if (b) update((s) => { s.layout = b.dataset.layout; });
  });
  $("clock-toggles").addEventListener("click", (e) => {
    const b = e.target.closest("[data-clock]");
    if (!b) return;
    update((s) => {
      const c = s.clocks.find((c) => c.id === Number(b.dataset.clock));
      if (c) c.visible = !c.visible;
    }, "analysis");
  });
  $("table-clock").addEventListener("change", (e) => update((s) => { s.tableClock = Number(e.target.value); }));
  $("table-dev").addEventListener("change", (e) => update((s) => { s.tableDev = e.target.value; }));

  function renderPicker() {
    const avail = ESTIMATORS.filter((e) => !state.devs.includes(e.id));
    const group = (axis, title) => {
      const items = avail.filter((e) => e.axis === axis);
      if (!items.length) return "";
      return `<div class="picker-group"><h3>${title}</h3>${items.map((e) =>
        `<button type="button" class="picker-item" data-add="${e.id}"><b>${e.label}</b><span>${escapeHtml(e.note)}</span></button>`).join("")}</div>`;
    };
    picker.innerHTML = (group("sigma_y", "Frequency σ<sub>y</sub>") + group("sigma_x", "Time σ<sub>x</sub>")) ||
      `<p class="muted small">Every estimator is already selected.</p>`;
  }

  // ── Refresh ──
  function refresh() {
    const colors = palette();
    const colorOf = (id) => colors[state.devs.indexOf(id) % colors.length];

    $("count").textContent = `${state.devs.length} selected`;
    $("chips").innerHTML = state.devs.length
      ? state.devs.map((id) => `<span class="chip" style="--chip:${colorOf(id)}">
          <span class="chip-swatch"></span>${label(id)}
          <button type="button" class="chip-x" data-remove="${id}" aria-label="Remove ${label(id)}">×</button></span>`).join("")
      : `<p class="muted small">No estimators yet. Add one below.</p>`;

    root.querySelector('[data-opt="grid"]').value = state.grid;
    root.querySelector('[data-opt="ci"]').checked = state.ci;
    const conf = root.querySelector('[data-opt="confidence"]');
    conf.value = String(state.confidence);
    conf.disabled = !state.ci;
    const band = root.querySelector('[data-opt="band"]');
    band.checked = state.band;
    band.disabled = !state.ci;
    root.querySelector('[data-opt="bias"]').checked = state.bias;

    const r = state.run;
    $("run").disabled = r.state === "running" || !state.clocks.length || !state.devs.length;
    $("last-run").textContent = r.lastMs !== null ? r.summary || "" : "";

    // Axis tabs with counts; fall back to the axis that has estimators.
    const nY = state.devs.filter((d) => byId[d]?.axis === "sigma_y").length;
    const nX = state.devs.length - nY;
    if (state.axis === "sigma_y" && nY === 0 && nX > 0) state.axis = "sigma_x";
    if (state.axis === "sigma_x" && nX === 0 && nY > 0) state.axis = "sigma_y";
    root.querySelectorAll("[data-axis]").forEach((b) => {
      b.classList.toggle("on", b.dataset.axis === state.axis);
      b.setAttribute("aria-checked", String(b.dataset.axis === state.axis));
      b.querySelector(".seg-n").textContent = b.dataset.axis === "sigma_y" ? (nY || "") : (nX || "");
    });
    const confText = state.ci ? `${num(state.confidence * 100, 3)} % confidence` : "no confidence intervals";
    $("plot-caption").textContent = `${confText} · ${GRIDS.find(([v]) => v === state.grid)[1].toLowerCase()} τ grid`;

    const multi = state.clocks.length > 1;
    $("layout-seg").hidden = !multi;
    $("clock-toggles").hidden = !multi;
    root.querySelectorAll("[data-layout]").forEach((b) => b.classList.toggle("on", b.dataset.layout === state.layout));
    $("clock-toggles").innerHTML = state.clocks.map((c, i) =>
      `<button type="button" class="clock-toggle ${c.visible ? "on" : ""}" data-clock="${c.id}" aria-pressed="${c.visible}">
        <svg width="22" height="8" aria-hidden="true"><line x1="1" y1="4" x2="21" y2="4" stroke="currentColor" stroke-width="2" stroke-dasharray="${dashArray(i)}"/></svg>${escapeHtml(c.name)}</button>`).join("");

    renderPlot(colorOf);
    renderTable();
  }

  function dashArray(i) {
    return ["", "6 3", "2 3", "6 3 2 3", "10 4"][i % 5];
  }

  function emptyMessage() {
    if (!state.core.ok) {
      return `<h3>No analysis core found</h3><p>Build <code>sigmatau-cli</code> (see <code>app/web/README.md</code>) or put <code>julia</code> on the PATH, then restart <code>serve.py</code>.</p>`;
    }
    if (!state.clocks.length) return `<h3>No data yet</h3><p>Add a clock on the <a href="#data">Data</a> page, or start with a synthetic record.</p><button class="btn btn-primary" type="button" data-action="synth">Use a synthetic record</button>`;
    if (!state.devs.length) return `<h3>No estimators selected</h3><p>Use <b>+ Add estimator</b> on the left.</p>`;
    const errs = activeClocks().filter((c) => c.error);
    if (errs.length) return `<h3>The analysis failed</h3>${errs.map((c) => `<p><b>${escapeHtml(c.name)}:</b> ${escapeHtml(c.error)}</p>`).join("")}`;
    if (!activeClocks().some((c) => c.result)) return state.run.state === "running" ? `<p class="muted">Running…</p>` : `<p class="muted">Press <b>Run analysis</b> to compute.</p>`;
    const other = state.axis === "sigma_y" ? "Time σx(τ)" : "Frequency σy(τ)";
    return `<p class="muted">None of the selected estimators plot on this axis. Try ${other}.</p>`;
  }

  function seriesFor(clock) {
    if (!clock.result) return [];
    return clock.result.results.filter((r) => r.axis === state.axis && state.devs.includes(r.dev));
  }

  function renderPlot(colorOf) {
    const Plotly = window.Plotly;
    const clocks = activeClocks();
    const side = state.layout === "side" && clocks.length > 1;
    const overlayMulti = !side && clocks.length > 1;
    const traces = [];
    clocks.forEach((clock, ci) => {
      const axisSuffix = side ? (ci === 0 ? "" : String(ci + 1)) : "";
      for (const r of seriesFor(clock)) {
        const color = colorOf(r.dev);
        const name = overlayMulti ? `${label(r.dev)} · ${clock.name}` : label(r.dev);
        const hasCI = state.ci && r.ci_lo.length === r.tau.length;
        const custom = r.tau.map((_, i) => [
          r.edf[i] ?? null, r.noise[i] ? noiseLabel(r.noise[i]) : "", r.ci_lo[i] ?? null, r.ci_hi[i] ?? null,
        ]);
        if (hasCI && state.band) {
          traces.push({
            x: [...r.tau, ...[...r.tau].reverse()],
            y: [...r.ci_hi, ...[...r.ci_lo].reverse()],
            fill: "toself", fillcolor: hexAlpha(color, 0.14), line: { width: 0 },
            hoverinfo: "skip", showlegend: false, type: "scatter", mode: "lines",
            xaxis: "x" + axisSuffix, yaxis: "y" + axisSuffix, legendgroup: name,
          });
        }
        traces.push({
          type: "scatter", mode: "lines+markers", name, legendgroup: name,
          showlegend: !side || ci === 0,
          x: r.tau, y: r.sigma, customdata: custom,
          xaxis: "x" + axisSuffix, yaxis: "y" + axisSuffix,
          line: { color, width: 2, dash: overlayMulti ? DASHES[ci % DASHES.length] : "solid" },
          marker: { color, size: 6, symbol: overlayMulti ? SYMBOLS[ci % SYMBOLS.length] : "circle" },
          error_y: hasCI && !state.band ? {
            type: "data", symmetric: false, thickness: 1.2, width: 3, color,
            array: r.ci_hi.map((h, i) => (h ?? r.sigma[i]) - r.sigma[i]),
            arrayminus: r.ci_lo.map((l, i) => r.sigma[i] - (l ?? r.sigma[i])),
          } : undefined,
          hovertemplate:
            `<b>${escapeHtml(name)}</b><br>τ = %{x:.4g} s<br>σ = %{y:.4e}` +
            (hasCI ? `<br>CI [%{customdata[2]:.3e}, %{customdata[3]:.3e}]<br>EDF %{customdata[0]:.4g} · %{customdata[1]}` : "") +
            "<extra></extra>",
        });
      }
    });

    const empty = $("empty");
    const showEmpty = traces.length === 0;
    empty.hidden = !showEmpty;
    plotEl.hidden = showEmpty;
    if (showEmpty) {
      empty.innerHTML = emptyMessage();
      if (Plotly) Plotly.purge(plotEl);
      return;
    }
    if (!Plotly) return;

    const layout = baseLayout();
    const xTitle = "Averaging time τ (s)";
    if (side) {
      layout.grid = { rows: 1, columns: clocks.length, pattern: "independent", xgap: 0.08 };
      layout.annotations = [];
      clocks.forEach((c, i) => {
        const k = i === 0 ? "" : String(i + 1);
        layout["xaxis" + k] = { ...layout._axis, ...logAxis(xTitle), ...(i ? { matches: "x" } : {}) };
        layout["yaxis" + k] = { ...layout._axis, ...logAxis(i === 0 ? AXIS_TITLE[state.axis] : ""),
          ...(i ? { matches: "y" } : {}) };
        layout.annotations.push({
          text: `<b>${escapeHtml(c.name)}</b>`, showarrow: false,
          xref: `x${k} domain`, yref: `y${k} domain`, x: 0.02, y: 0.02, xanchor: "left", yanchor: "bottom",
          font: { size: 13 },
        });
      });
    } else {
      layout.xaxis = { ...layout._axis, ...logAxis(xTitle) };
      layout.yaxis = { ...layout._axis, ...logAxis(AXIS_TITLE[state.axis]) };
    }
    delete layout._axis;
    Plotly.react(plotEl, traces, layout, CONFIG);
  }

  function renderTable() {
    const withResult = state.clocks.filter((c) => c.result);
    const clockSel = $("table-clock"), devSel = $("table-dev");
    if (!withResult.some((c) => c.id === state.tableClock)) state.tableClock = withResult[0]?.id ?? null;
    clockSel.hidden = withResult.length < 2;
    clockSel.innerHTML = withResult.map((c) => `<option value="${c.id}">${escapeHtml(c.name)}</option>`).join("");
    if (state.tableClock !== null) clockSel.value = String(state.tableClock);

    const clock = withResult.find((c) => c.id === state.tableClock);
    const results = clock ? clock.result.results.filter((r) => state.devs.includes(r.dev)) : [];
    if (!results.some((r) => r.dev === state.tableDev)) {
      state.tableDev = (results.find((r) => r.axis === state.axis) ?? results[0])?.dev ?? null;
    }
    devSel.innerHTML = results.map((r) => `<option value="${r.dev}">${label(r.dev)}</option>`).join("");
    devSel.hidden = !results.length;
    if (state.tableDev) devSel.value = state.tableDev;

    const r = results.find((x) => x.dev === state.tableDev);
    const table = $("table");
    if (!r) {
      table.innerHTML = "";
      $("table-foot").textContent = "Results appear here after a run.";
      return;
    }
    const hasCI = r.ci_lo.length === r.tau.length;
    // The interval bounds go in each row's tooltip and in the TSV export.
    table.innerHTML = `<thead><tr><th>τ (s)</th><th>σ</th>${hasCI ? "<th>EDF</th><th>noise</th>" : ""}</tr></thead>
      <tbody>${r.tau.map((t, i) => `<tr${hasCI ? ` title="CI ${sci(r.ci_lo[i], 3)} – ${sci(r.ci_hi[i], 3)}"` : ""}><td>${num(t, 5)}</td><td>${sci(r.sigma[i])}</td>${hasCI
        ? `<td>${num(r.edf[i], 4)}</td><td><span class="badge">${noiseLabel(r.noise[i])}</span></td>`
        : ""}</tr>`).join("")}</tbody>`;
    const rec = clock.result.record;
    $("table-foot").textContent = `${label(r.dev)} on ${clock.name}: ${r.tau.length} τ values from ${thousands(rec.n)} samples.`;
  }

  return { refresh };
}

function hexAlpha(hex, a) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}
