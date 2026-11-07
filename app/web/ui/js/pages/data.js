// Data page: the list of clocks (records) and the import dialog. Each clock
// keeps its own reader settings; the core reads the file with them.
import { state, update, newClock } from "../state.js";
import { sniff, tau0FromColumn } from "../import.js";
import { baseLayout, CONFIG, palette, DASHES } from "../plot.js";
import { escapeHtml, thousands, duration, num } from "../format.js";
import { openFile } from "../platform.js";

const SCALES = [[1, "seconds (s)"], [1e-9, "nanoseconds (ns)"], [1e-12, "picoseconds (ps)"]];
const DETRENDS = [["none", "None"], ["mean", "Remove mean"], ["endpoint", "Endpoint line"], ["linear", "Linear fit"]];

export function mountData(root) {
  root.innerHTML = `
  <div class="data-page">
    <section class="clocks">
      <div class="section-head"><h2>Clocks</h2>
        <button class="btn btn-primary" type="button" data-ref="add">+ Add clock</button></div>
      <div class="clock-list" data-ref="list"></div>
    </section>
    <section class="panel record-card">
      <div class="section-head"><h2>Record</h2><span class="muted small" data-ref="record-caption"></span></div>
      <div class="plot record-plot" data-ref="plot"></div>
      <div class="empty" data-ref="empty" hidden><p class="muted">The phase record appears here after a run.</p></div>
    </section>
  </div>
  <dialog class="dialog import-dialog" data-ref="dialog"></dialog>`;

  const $ = (r) => root.querySelector(`[data-ref="${r}"]`);
  const dialog = $("dialog");

  $("add").addEventListener("click", () => openImport(dialog));
  $("list").addEventListener("click", (e) => {
    const b = e.target.closest("[data-act]");
    if (!b) return;
    const id = Number(b.closest("[data-id]").dataset.id);
    const clock = state.clocks.find((c) => c.id === id);
    if (b.dataset.act === "remove") update((s) => { s.clocks = s.clocks.filter((c) => c.id !== id); }, "analysis");
    if (b.dataset.act === "edit") openImport(dialog, clock);
  });
  $("list").addEventListener("change", (e) => {
    const inp = e.target.closest("[data-rename]");
    if (!inp) return;
    const id = Number(inp.closest("[data-id]").dataset.id);
    update((s) => { s.clocks.find((c) => c.id === id).name = inp.value.trim() || "Clock"; });
  });

  function refresh() {
    const colors = palette();
    $("list").innerHTML = state.clocks.length ? state.clocks.map((c, i) => {
      const rec = c.result?.record;
      const src = c.source.type === "synth" ? `Synthetic, ${thousands(c.source.n)} samples, seed ${c.source.seed}` : escapeHtml(c.source.name);
      const facts = rec
        ? `${thousands(rec.n)} samples · τ₀ ${num(rec.tau0, 4)} s · ${duration(rec.n * rec.tau0)} · ${c.kind} · detrend ${c.detrend}`
        : `${c.kind} · detrend ${c.detrend}${c.error ? "" : " · not run yet"}`;
      return `<article class="clock-card panel" data-id="${c.id}" style="--chip:${colors[i % colors.length]}">
        <input class="clock-name" data-rename value="${escapeHtml(c.name)}" aria-label="Clock name">
        <p class="small">${src}</p>
        <p class="muted small">${facts}</p>
        ${c.error ? `<p class="error small">${escapeHtml(c.error)}</p>` : ""}
        <div class="card-actions">
          <button class="btn btn-ghost" type="button" data-act="edit">Edit import</button>
          <button class="btn btn-ghost" type="button" data-act="remove">Remove</button>
        </div></article>`;
    }).join("") : `<div class="empty-card panel"><p>No clocks yet. Import a phase or frequency file, or start with a synthetic record.</p></div>`;

    const withSeries = state.clocks.filter((c) => c.visible && c.result?.record?.series?.t?.length);
    const Plotly = window.Plotly;
    $("empty").hidden = withSeries.length > 0;
    $("plot").hidden = withSeries.length === 0;
    $("record-caption").textContent = withSeries.length ? "Phase x(t), decimated for display" : "";
    if (!Plotly || !withSeries.length) return;
    const traces = withSeries.map((c, i) => ({
      type: "scattergl", mode: "lines", name: c.name,
      x: c.result.record.series.t, y: c.result.record.series.x,
      line: { color: colors[state.clocks.indexOf(c) % colors.length], width: 1.2, dash: DASHES[0] },
      hovertemplate: `<b>${escapeHtml(c.name)}</b><br>t = %{x:.4g} s<br>x = %{y:.4e} s<extra></extra>`,
    }));
    const layout = baseLayout();
    layout.xaxis = { ...layout._axis, title: { text: "Time t (s)" } };
    layout.yaxis = { ...layout._axis, title: { text: "Phase x (s)" }, exponentformat: "power" };
    delete layout._axis;
    layout.showlegend = withSeries.length > 1;
    Plotly.react($("plot"), traces, layout, CONFIG);
  }

  return { refresh };
}

// ── Import dialog ─────────────────────────────────────────────────────────

function openImport(dialog, existing = null) {
  const draft = existing
    ? { ...existing, source: { ...existing.source } }
    : newClock({ source: { type: "file", name: "", content: "" }, value_col: 1, time_col: 0 });
  let preview = draft.source.type === "file" && draft.source.content ? sniff(draft.source.name, draft.source.content) : null;
  let mode = draft.source.type;

  function render() {
    const colCount = preview?.nCols ?? 0;
    const role = (k) => (k === draft.value_col ? "value" : k === draft.time_col ? "time" : "");
    const head = Array.from({ length: colCount }, (_, i) => {
      const k = i + 1, r = role(k);
      const nm = preview.headerNames[i] ? escapeHtml(preview.headerNames[i]) : `Column ${k}`;
      return `<th><button type="button" class="col-head ${r}" data-col="${k}" title="Click to cycle: value, time, unused">
        <span>${nm}</span><span class="role">${r === "value" ? "Value" : r === "time" ? "Time" : "unused"}</span></button></th>`;
    }).join("");
    const body = (preview?.rows ?? []).map((row) =>
      `<tr>${Array.from({ length: colCount }, (_, i) => `<td class="${role(i + 1)}">${escapeHtml(row[i] ?? "")}</td>`).join("")}</tr>`).join("");
    const tauHint = draft.time_col > 0 && preview ? tau0FromColumn(preview.rows, draft.time_col) : preview?.tau0Hint;

    dialog.innerHTML = `<form method="dialog" class="import">
      <h2>${existing ? "Edit import" : "Add a clock"}</h2>
      <div class="seg" role="radiogroup">
        <button type="button" data-mode="file" class="${mode === "file" ? "on" : ""}">From a file</button>
        <button type="button" data-mode="synth" class="${mode === "synth" ? "on" : ""}">Synthetic record</button>
      </div>
      ${mode === "file" ? `
        <button type="button" class="drop" data-ref="drop">
          <span>${draft.source.name ? `<b>${escapeHtml(draft.source.name)}</b> · ${thousands(preview?.nRows ?? 0)} rows · ${colCount} columns` : "Drop a file here or <u>choose one</u>"}</span>
        </button>
        ${preview ? `<p class="muted small">Click a column header to set it as the <b>Value</b> column, click again for <b>Time</b>, again to leave it unused.
          ${preview.comments ? `${preview.comments} comment line(s) and ` : ""}${draft.header} header line(s) skipped.</p>
          <div class="preview-wrap"><table class="preview"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>` : ""}
      ` : `
        <div class="grid2">
          <label class="field"><span>Samples</span><input type="number" min="16" step="1" data-f="synth_n" value="${draft.source.n ?? 16384}"></label>
          <label class="field"><span>Seed</span><input type="number" step="1" data-f="synth_seed" value="${draft.source.seed ?? 42}"></label>
        </div>
        <p class="muted small">White FM plus random-walk FM, for trying the explorer without a file.</p>`}
      <div class="grid2">
        <label class="field"><span>Name</span><input type="text" data-f="name" value="${escapeHtml(draft.name)}"></label>
        ${mode === "file" ? `<label class="field"><span>Values are</span><select data-f="kind">
          <option value="phase" ${draft.kind === "phase" ? "selected" : ""}>Phase x(t)</option>
          <option value="frequency" ${draft.kind === "frequency" ? "selected" : ""}>Fractional frequency y(t)</option></select></label>` : "<span></span>"}
        <label class="field"><span>τ₀ (s)</span><input type="number" step="any" min="0" data-f="tau0" value="${draft.tau0 ?? ""}" placeholder="${tauHint ? `auto: ${num(tauHint, 4)}` : mode === "synth" ? "1" : "from file"}"></label>
        ${mode === "file" ? `<label class="field"><span>Units</span><select data-f="scaling">${SCALES.map(([v, t]) =>
          `<option value="${v}" ${draft.scaling === v ? "selected" : ""}>${t}</option>`).join("")}</select></label>` : "<span></span>"}
        <label class="field"><span>Detrend</span><select data-f="detrend">${DETRENDS.map(([v, t]) =>
          `<option value="${v}" ${draft.detrend === v ? "selected" : ""}>${t}</option>`).join("")}</select></label>
      </div>
      <p class="error small" data-ref="err"></p>
      <div class="dialog-actions">
        <button class="btn" value="cancel">Cancel</button>
        <button class="btn btn-primary" type="button" data-ref="ok">${existing ? "Save" : "Add clock"}</button>
      </div></form>`;
  }

  function takeFile(file) {
    if (!file) return;
    const { name, text } = file;
    draft.source = { type: "file", name, content: text };
    preview = sniff(name, text);
    draft.header = preview.header;
    // Default roles: one column = value; two or more = time then value.
    if (preview.nCols >= 2) { draft.time_col = 1; draft.value_col = 2; } else { draft.time_col = 0; draft.value_col = 1; }
    if (!existing || draft.name.startsWith("Clock ")) draft.name = name.replace(/\.[^.]+$/, "");
    render();
  }

  dialog.onclick = (e) => {
    const t = e.target;
    const modeBtn = t.closest("[data-mode]");
    if (modeBtn) {
      mode = modeBtn.dataset.mode;
      if (mode === "synth") draft.source = { type: "synth", n: draft.source.n ?? 16384, seed: draft.source.seed ?? 42 };
      else if (draft.source.type !== "file") { draft.source = { type: "file", name: "", content: "" }; preview = null; }
      return render();
    }
    if (t.closest('[data-ref="drop"]')) {
      openFile().then(takeFile);
      return;
    }
    const col = t.closest("[data-col]");
    if (col) {
      const k = Number(col.dataset.col);
      if (draft.value_col === k) { draft.value_col = 0; draft.time_col = k; }
      else if (draft.time_col === k) { draft.time_col = 0; }
      else draft.value_col = k;   // the previous value column becomes unused
      return render();
    }
    if (t.closest('[data-ref="ok"]')) {
      const err = dialog.querySelector('[data-ref="err"]');
      if (mode === "file" && !draft.source.content) { err.textContent = "Choose a file first."; return; }
      if (mode === "file" && !draft.value_col) { err.textContent = "Pick a Value column by clicking its header."; return; }
      if (mode === "file" && !draft.time_col && !draft.tau0 && !preview?.tau0Hint) { err.textContent = "Without a time column, enter τ₀."; return; }
      update((s) => {
        if (existing) Object.assign(s.clocks.find((c) => c.id === existing.id), draft, { result: null, error: null });
        else s.clocks.push(draft);
      }, "analysis");
      dialog.close();
    }
  };
  dialog.onchange = (e) => {
    const t = e.target;
    const f = t.dataset.f;
    if (!f) return;
    if (f === "synth_n") draft.source.n = Math.max(16, Math.round(Number(t.value)));
    else if (f === "synth_seed") draft.source.seed = Math.round(Number(t.value));
    else if (f === "tau0") draft.tau0 = t.value === "" ? null : Number(t.value);
    else if (f === "scaling") draft.scaling = Number(t.value);
    else draft[f] = t.value;
  };
  dialog.ondragover = (e) => { e.preventDefault(); dialog.querySelector(".drop")?.classList.add("over"); };
  dialog.ondragleave = () => dialog.querySelector(".drop")?.classList.remove("over");
  dialog.ondrop = async (e) => {
    e.preventDefault();
    const f = e.dataTransfer.files[0];
    if (f) takeFile({ name: f.name, text: await f.text() });
  };

  render();
  dialog.showModal();
}
