// The import dialog: shows the top of a file, lets the user say which column
// is time and which is the value, the record type, units and τ0, and returns
// the reader settings for the core. openImport() resolves to null on cancel.
import { sniff, guessRoles, timeStats, thin, PHASE_UNITS, FREQ_UNITS, TIME_UNITS } from "./import.js";
import { escapeHtml as h, thousands, duration, num } from "./format.js";
import { seriesColor, hexAlpha } from "./plot.js";

const FILE_ICON = `<svg width="22" height="26" viewBox="0 0 22 26" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><path d="M3 2h10l6 6v16H3z"/><path d="M13 2v6h6"/></svg>`;

// Settings are remembered per file format: same extension, column count
// and column names.
const formatKey = (s) => `sigmatau.import.${s.name.toLowerCase().split(".").pop()}.${s.nCols}.${s.names.join("|")}`;
function remembered(s) {
  try { return JSON.parse(localStorage.getItem(formatKey(s)) ?? "null"); } catch { return null; }
}
function remember(s, st) {
  try {
    if (!st.remember) { localStorage.removeItem(formatKey(s)); return; }
    localStorage.setItem(formatKey(s), JSON.stringify({ kind: st.kind, time_col: st.time, timeUnit: st.timeUnit, value_cols: st.values, scaling: st.scaling, tau0: st.tau0In, against: st.against }));
  } catch { /* storage off: nothing to remember */ }
}

// `prev`: settings of the clocks being re-imported (Import settings…); it
// always shows the dialog. Otherwise settings remembered for files like this
// one load the file without asking.
export function openImport(file, prev = null) {
  const dlg = document.getElementById("import-dialog");
  const s = sniff(file.name, file.text);
  const guess = guessRoles(s);
  const mem = prev ? null : remembered(s);
  const p = prev ?? mem;
  const st = {
    kind: p?.kind ?? "phase",
    time: p?.time_col ?? guess.time,          // 0 = no time column
    timeUnit: p?.timeUnit ?? guess.timeUnit,
    values: (p?.value_cols ?? [guess.value]).filter((c) => c >= 1 && c <= Math.max(1, s.nCols)),
    scaling: p?.scaling ?? 1,
    tau0In: p?.tau0 ?? s.tau0Hint ?? "",     // no guess: without a time column τ0 must be given
    against: p?.against ?? "",
    remember: true,
  };
  if (!st.values.length) st.values = [guess.value];
  st.values = st.values.filter((c) => c !== st.time);
  if (!st.values.length) st.time = 0, st.values = [guess.value];
  const fromMemory = !!mem;

  const derived = () => {
    const ts = st.time ? timeStats(s.columns[st.time - 1] ?? [], st.timeUnit) : null;
    const tau0 = ts && ts.tau0 > 0 ? ts.tau0 : Number(st.tau0In);
    return { ts, tau0 };
  };
  const colName = (c) => s.names[c - 1] || `column ${c}`;
  const settings = () => {
    const { ts, tau0 } = derived();
    const units = st.kind === "phase" ? PHASE_UNITS : FREQ_UNITS;
    return {
      kind: st.kind,
      time_col: st.time,
      timeUnit: st.timeUnit,
      value_cols: [...st.values],
      names: st.values.map((c) => s.names[c - 1] || null),
      against: st.against.trim(),
      header: s.header,
      scaling: st.scaling,
      tau0,
      hasTimeCol: st.time > 0,
      gaps: ts ? ts.gaps : 0,
      missing: ts ? ts.missing : 0,         // samples missing from the time column
      sourceInfo: `${s.stable32 ? "Stable32 file" : s.delimName.replace("-separated", "")} · ${st.kind} · τ₀ ${num(tau0, 4)} s` +
        (s.comments + s.header ? ` · ${s.comments + s.header} header line${s.comments + s.header > 1 ? "s" : ""} read automatically` : ""),
      unitName: (units.find(([f]) => f === st.scaling) ?? units[0])[1],
      remembered: fromMemory && !prev,
    };
  };

  // Remembered format: load straight away when the file still makes sense.
  if (fromMemory && !prev && s.nCols && !s.bad && derived().tau0 > 0) return Promise.resolve(settings());

  return new Promise((resolve) => {
    const finish = (out) => { dlg.close(); resolve(out); };

    const render = () => {
      const { ts, tau0 } = derived();
      const units = st.kind === "phase" ? PHASE_UNITS : FREQ_UNITS;
      if (!units.some(([f]) => f === st.scaling)) st.scaling = 1;
      const unitName = units.find(([f]) => f === st.scaling)[1];
      const role = (c) => (c === st.time ? "time" : st.values.includes(c) ? "value" : "ignore");
      const several = st.values.length > 1;
      const color = (c) => seriesColor(st.values.indexOf(c));

      const pills = Array.from({ length: s.nCols }, (_, i) => {
        const c = i + 1, r = role(c);
        const timeOpts = TIME_UNITS.map(([f, l]) =>
          `<option value="time:${f}" ${r === "time" && f === st.timeUnit ? "selected" : ""}>${l} ▾</option>`).join("");
        return `<select class="${r}" data-col="${c}" aria-label="Column ${c} holds" ${r === "value" ? `style="background:${color(c)};border-color:${color(c)}"` : ""}>
          <option value="ignore" ${r === "ignore" ? "selected" : ""}>Ignore ▾</option>
          ${timeOpts}
          <option value="value" ${r === "value" ? "selected" : ""}>Clock: ${h(colName(c))} ▾</option>
        </select>`;
      }).join("");

      const cls = (c) => (role(c) === "ignore" ? "" : role(c));
      const tint = (c, a) => (role(c) === "value" ? ` style="background:${hexAlpha(color(c), a)}${a > 0.15 ? `;color:${color(c)}` : ""}"` : "");
      const topRows = s.top.map(({ ln, text }) => {
        if (ln === s.namesLine) {
          return `<tr class="names"><td class="ln">${ln}</td>${s.names.map((n, i) =>
            `<td class="${cls(i + 1)}"${tint(i + 1, 0.2)} data-col="${i + 1}" title="Click to add or drop this column as a clock">${h(n)}</td>`).join("")}</tr>`;
        }
        return `<tr class="comment"><td class="ln">${ln}</td><td colspan="${Math.max(1, s.nCols)}">${h(text)}</td></tr>`;
      }).join("");
      const bodyRows = s.rows.map(({ ln, f }) =>
        `<tr><td class="ln">${ln}</td>${Array.from({ length: s.nCols }, (_, i) =>
          `<td class="${cls(i + 1)}"${tint(i + 1, 0.1)}>${h(f[i] ?? "")}</td>`).join("")}</tr>`).join("");
      const more = s.nRows - s.rows.length;

      const n = s.nRows + (ts ? ts.missing : 0);
      const spanS = ts ? ts.span : tau0 * (s.nRows - 1);
      const summary = s.nCols
        ? `${several ? `${st.values.length} clocks · ${thousands(s.nRows)} samples each` : `${thousands(s.nRows)} samples`} · ${duration(spanS)} · ${st.kind === "phase" ? "phase" : "frequency"} in ${unitName}${st.against.trim() ? ` · vs ${h(st.against.trim())}` : ""}`
        : "";
      const err = !s.nCols ? "No numeric rows found in this file."
        : s.bad ? `Line ${s.bad} is not numeric. The core will refuse the file.`
        : !(tau0 > 0) ? "Enter the sample interval τ₀: this file has no time column." : "";

      // One normalised line per clock column.
      const lines = st.values.map((c, k) => {
        const pts = thin(s.columns[c - 1] ?? [], 200).filter((p) => Number.isFinite(p[1]));
        if (pts.length < 2) return "";
        const ys = pts.map((p) => p[1]);
        const lo = Math.min(...ys), hi = Math.max(...ys), x1 = pts[pts.length - 1][0] || 1;
        const sy = (y) => (hi > lo ? 80 - ((y - lo) / (hi - lo)) * 64 : 48);
        return `<polyline fill="none" stroke="${seriesColor(k)}" stroke-width="1.6" vector-effect="non-scaling-stroke" points="${pts.map((p) => `${(12 + (p[0] / x1) * 306).toFixed(1)},${sy(p[1]).toFixed(1)}`).join(" ")}"/>`;
      }).join("");
      const spark = lines ? `<svg viewBox="0 0 330 92" preserveAspectRatio="none" width="100%" height="100%">${lines}</svg>` : "";

      const gap = ts && ts.gaps
        ? `<div class="gap-card"><b>${ts.gaps} gap${ts.gaps > 1 ? "s" : ""} in the time column</b>${thousands(ts.missing)} missing sample${ts.missing === 1 ? "" : "s"}. Turn on <i>Fill gaps</i> on the Data page to put them back on an even grid.</div>`
        : "";
      const tauLine = ts && ts.tau0 > 0
        ? `<div class="tau-line"><span class="mono">${num(ts.tau0, 4)} s</span><span class="ok">✓ from the time column${ts.even ? ", evenly spaced" : ", median step"}</span></div>`
        : `<div class="tau-line"><input type="number" id="imp-tau0" min="0" step="any" value="${h(st.tau0In)}"> <span class="muted small">s${s.tau0Hint ? " · from the file header" : " · no time column"}</span></div>`;

      dlg.innerHTML = `
        <div class="imp-head">${FILE_ICON}
          <div><b>Import ${h(s.name)}</b><span class="muted small">${thousands(s.lines)} lines · ${s.delimName}${s.stable32 ? " · Stable32 file" : ""}</span></div>
          <span class="grow"></span><span>Click a column heading to choose what it holds</span>
        </div>
        <div class="imp-body">
          <div class="imp-main">
            <div class="roles" style="grid-template-columns: repeat(${Math.max(1, s.nCols)}, minmax(0, 1fr))">${pills}</div>
            <div class="preview"><table>
              <colgroup><col style="width:44px">${"<col>".repeat(Math.max(1, s.nCols))}</colgroup>
              ${topRows}${bodyRows}
              ${more > 0 ? `<tr><td class="more" colspan="${s.nCols + 1}">… ${thousands(more)} more rows</td></tr>` : ""}
            </table></div>
            <div class="imp-notes">
              ${s.comments ? `<span>Skipped: lines starting with #</span>` : ""}
              ${s.namesLine ? `<span>Column names from line ${s.namesLine}</span>` : s.header ? `<span>${s.header} header line${s.header > 1 ? "s" : ""} skipped</span>` : ""}
            </div>
          </div>
          <div class="imp-side">
            <div><label class="lbl">Record type</label>
              <div class="seg" role="group"><button type="button" data-kind="phase" class="${st.kind === "phase" ? "on" : ""}">Phase (time error)</button><button type="button" data-kind="frequency" class="${st.kind === "frequency" ? "on" : ""}">Fractional frequency</button></div></div>
            <div><label class="lbl" for="imp-units">${st.kind === "phase" ? "Phase units" : "Frequency units"}</label>
              <select id="imp-units">${units.map(([f, l]) => `<option value="${f}" ${f === st.scaling ? "selected" : ""}>${l}</option>`).join("")}</select></div>
            <div><label class="lbl">Sample interval τ<sub>0</sub></label>${tauLine}</div>
            ${gap}
            <div><label class="lbl" for="imp-against">Measured against</label>
              <input type="text" id="imp-against" class="wide-input" placeholder="e.g. H-maser" value="${h(st.against)}">
              <div class="muted small" style="margin-top:6px">You can switch the reference to any clock later, in the Clocks list.</div></div>
            <div><label class="lbl">${several ? `${st.values.length} clocks from this file` : `Preview of ${h(colName(st.values[0]))}`}</label><div class="mini-plot">${spark}</div></div>
            <label class="check"><input type="checkbox" id="imp-remember" ${st.remember ? "checked" : ""}>Remember for files like this</label>
          </div>
        </div>
        <div class="imp-foot">
          <span class="summary">${err ? `<span class="err">${h(err)}</span>` : h(summary)}</span>
          <button class="btn" type="button" data-act="cancel">Cancel</button>
          <button class="btn btn-primary" type="button" data-act="load" ${err ? "disabled" : ""}>${several ? `Load ${st.values.length} clocks` : "Load record"}</button>
        </div>`;
    };

    dlg.onclick = (e) => {
      const kind = e.target.closest("[data-kind]")?.dataset.kind;
      if (kind) { st.kind = kind; render(); return; }
      const col = e.target.closest("tr.names td[data-col]")?.dataset.col;
      if (col) {
        const c = Number(col);
        if (st.values.includes(c)) { if (st.values.length > 1) st.values = st.values.filter((v) => v !== c); }
        else { st.values.push(c); st.values.sort((a, b) => a - b); if (st.time === c) st.time = 0; }
        render(); return;
      }
      const act = e.target.closest("[data-act]")?.dataset.act;
      if (act === "cancel") finish(null);
      if (act === "load") {
        remember(s, st);
        finish(settings());
      }
    };
    dlg.onchange = (e) => {
      const t = e.target;
      if (t.dataset.col) {
        const c = Number(t.dataset.col), v = t.value;
        if (v === "value") { if (!st.values.includes(c)) st.values.push(c); if (st.time === c) st.time = 0; }
        else if (v.startsWith("time:")) {
          st.time = c; st.timeUnit = Number(v.slice(5));
          st.values = st.values.filter((x) => x !== c);
          if (!st.values.length) st.values = [Array.from({ length: s.nCols }, (_, i) => i + 1).find((x) => x !== c) ?? 1];
        } else {
          if (st.time === c) st.time = 0;
          if (st.values.includes(c) && st.values.length > 1) st.values = st.values.filter((x) => x !== c);
        }
        st.values.sort((a, b) => a - b);
        render();
      } else if (t.id === "imp-units") { st.scaling = Number(t.value); render(); }
      else if (t.id === "imp-remember") st.remember = t.checked;
    };
    dlg.oninput = (e) => {
      if (e.target.id === "imp-against") { st.against = e.target.value; return; }
      if (e.target.id === "imp-tau0") {
        st.tau0In = e.target.value;
        const { tau0 } = derived();
        dlg.querySelector("[data-act=load]").disabled = !(tau0 > 0);
        if (tau0 > 0) dlg.querySelector(".imp-foot .err")?.remove();
      }
    };
    dlg.oncancel = (e) => { e.preventDefault(); finish(null); };

    render();
    dlg.showModal();
  });
}
