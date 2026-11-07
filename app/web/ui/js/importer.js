// The import dialog: shows the top of a file, lets the user say which column
// is time and which is the value, the record type, units and τ0, and returns
// the reader settings for the core. openImport() resolves to null on cancel.
import { sniff, guessRoles, timeStats, thin, PHASE_UNITS, FREQ_UNITS, TIME_UNITS } from "./import.js";
import { escapeHtml as h, thousands, duration, num } from "./format.js";

const FILE_ICON = `<svg width="22" height="26" viewBox="0 0 22 26" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><path d="M3 2h10l6 6v16H3z"/><path d="M13 2v6h6"/></svg>`;

export function openImport(file, prev = null) {
  const dlg = document.getElementById("import-dialog");
  const s = sniff(file.name, file.text);
  const guess = guessRoles(s);
  const st = {
    kind: prev?.kind ?? "phase",
    time: prev?.time_col ?? guess.time,          // 0 = no time column
    timeUnit: prev?.timeUnit ?? guess.timeUnit,
    value: Math.min(prev?.value_col ?? guess.value, Math.max(1, s.nCols)),
    scaling: prev?.scaling ?? 1,
    tau0In: prev?.tau0 ?? s.tau0Hint ?? 1,
  };
  if (st.time === st.value) st.time = 0;

  return new Promise((resolve) => {
    const finish = (out) => { dlg.close(); resolve(out); };

    const derived = () => {
      const ts = st.time ? timeStats(s.columns[st.time - 1] ?? [], st.timeUnit) : null;
      const tau0 = ts && ts.tau0 > 0 ? ts.tau0 : Number(st.tau0In);
      return { ts, tau0 };
    };

    const render = () => {
      const { ts, tau0 } = derived();
      const units = st.kind === "phase" ? PHASE_UNITS : FREQ_UNITS;
      if (!units.some(([f]) => f === st.scaling)) st.scaling = 1;
      const unitName = units.find(([f]) => f === st.scaling)[1];
      const valueName = s.names[st.value - 1] || `column ${st.value}`;
      const role = (c) => (c === st.time ? "time" : c === st.value ? "value" : "ignore");

      const pills = Array.from({ length: s.nCols }, (_, i) => {
        const c = i + 1, r = role(c);
        const timeOpts = TIME_UNITS.map(([f, l]) =>
          `<option value="time:${f}" ${r === "time" && f === st.timeUnit ? "selected" : ""}>${l} ▾</option>`).join("");
        return `<select class="${r}" data-col="${c}" aria-label="Column ${c} holds">
          <option value="ignore" ${r === "ignore" ? "selected" : ""}>Ignore ▾</option>
          ${timeOpts}
          <option value="value" ${r === "value" ? "selected" : ""}>Value: ${st.kind === "phase" ? "phase" : "frequency"} ▾</option>
        </select>`;
      }).join("");

      const cls = (c) => (role(c) === "ignore" ? "" : role(c));
      const topRows = s.top.map(({ ln, text }) => {
        if (ln === s.namesLine) {
          return `<tr class="names"><td class="ln">${ln}</td>${s.names.map((n, i) =>
            `<td class="${cls(i + 1)}" data-col="${i + 1}" title="Use as the value column">${h(n)}</td>`).join("")}</tr>`;
        }
        return `<tr class="comment"><td class="ln">${ln}</td><td colspan="${Math.max(1, s.nCols)}">${h(text)}</td></tr>`;
      }).join("");
      const bodyRows = s.rows.map(({ ln, f }) =>
        `<tr><td class="ln">${ln}</td>${Array.from({ length: s.nCols }, (_, i) =>
          `<td class="${cls(i + 1)}">${h(f[i] ?? "")}</td>`).join("")}</tr>`).join("");
      const more = s.nRows - s.rows.length;

      const n = s.nRows + (ts ? ts.missing : 0);
      const spanS = ts ? ts.span : tau0 * (s.nRows - 1);
      const summary = s.nCols
        ? `${thousands(s.nRows)} samples · ${duration(spanS)} · ${st.kind === "phase" ? "phase" : "frequency"} in ${unitName}`
        : "";
      const err = !s.nCols ? "No numeric rows found in this file."
        : s.bad ? `Line ${s.bad} is not numeric. The core will refuse the file.`
        : !(tau0 > 0) ? "Enter a sample interval τ₀ greater than zero." : "";

      const pts = thin(s.columns[st.value - 1] ?? []);
      let spark = "";
      if (pts.length > 1) {
        const ys = pts.map((p) => p[1]).filter(Number.isFinite);
        const lo = Math.min(...ys), hi = Math.max(...ys), x1 = pts[pts.length - 1][0] || 1;
        const sy = (y) => (hi > lo ? 80 - ((y - lo) / (hi - lo)) * 64 : 48);
        spark = `<svg viewBox="0 0 330 92" preserveAspectRatio="none" width="100%" height="100%"><polyline fill="none" stroke="var(--s1)" stroke-width="1.6" vector-effect="non-scaling-stroke" points="${pts.filter((p) => Number.isFinite(p[1])).map((p) => `${(12 + (p[0] / x1) * 306).toFixed(1)},${sy(p[1]).toFixed(1)}`).join(" ")}"/></svg>`;
      }

      const gap = ts && ts.gaps
        ? `<div class="gap-card"><b>${ts.gaps} gap${ts.gaps > 1 ? "s" : ""} in the time column</b>${thousands(ts.missing)} missing sample${ts.missing === 1 ? "" : "s"}. Turn on <i>Fill gaps</i> on the Data page to put them back on an even grid.</div>`
        : "";
      const tauLine = ts && ts.tau0 > 0
        ? `<div class="tau-line"><span class="mono">${num(ts.tau0, 4)} s</span><span class="ok">✓ from the time column${ts.even ? ", evenly spaced" : ", median step"}</span></div>`
        : `<div class="tau-line"><input type="number" id="imp-tau0" min="0" step="any" value="${h(st.tau0In)}"> <span class="muted small">s${s.tau0Hint ? " · from the file header" : " · no time column"}</span></div>`;

      dlg.innerHTML = `
        <div class="imp-head">${FILE_ICON}
          <div><b>Import ${h(s.name)}</b><span class="muted small">${thousands(s.lines)} lines · ${s.delimName}${s.stable32 ? " · Stable32 file" : ""}</span></div>
          <span class="grow"></span><span>Choose what each column holds with the menus above it</span>
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
            <div><label class="lbl">Preview of ${h(valueName)}</label><div class="mini-plot">${spark}</div></div>
          </div>
        </div>
        <div class="imp-foot">
          <span class="summary">${err ? `<span class="err">${h(err)}</span>` : h(summary)}</span>
          <button class="btn" type="button" data-act="cancel">Cancel</button>
          <button class="btn btn-primary" type="button" data-act="load" ${err ? "disabled" : ""}>Load record</button>
        </div>`;
    };

    dlg.onclick = (e) => {
      const kind = e.target.closest("[data-kind]")?.dataset.kind;
      if (kind) { st.kind = kind; render(); return; }
      const col = e.target.closest("tr.names td[data-col]")?.dataset.col;
      if (col) { st.value = Number(col); if (st.time === st.value) st.time = 0; render(); return; }
      const act = e.target.closest("[data-act]")?.dataset.act;
      if (act === "cancel") finish(null);
      if (act === "load") {
        const { ts, tau0 } = derived();
        const units = st.kind === "phase" ? PHASE_UNITS : FREQ_UNITS;
        finish({
          kind: st.kind,
          time_col: st.time,
          timeUnit: st.timeUnit,
          value_col: st.value,
          header: s.header,
          scaling: st.scaling,
          tau0,
          hasTimeCol: st.time > 0,
          gaps: ts ? ts.gaps : 0,
          sourceInfo: `${s.stable32 ? "Stable32 file" : s.delimName.replace("-separated", "")} · ${st.kind} · τ₀ ${num(tau0, 4)} s` +
            (s.comments + s.header ? ` · ${s.comments + s.header} header line${s.comments + s.header > 1 ? "s" : ""} read automatically` : ""),
          unitName: units.find(([f]) => f === st.scaling)[1],
        });
      }
    };
    dlg.onchange = (e) => {
      const t = e.target;
      if (t.dataset.col) {
        const c = Number(t.dataset.col), v = t.value;
        if (v === "value") { st.value = c; if (st.time === c) st.time = 0; }
        else if (v.startsWith("time:")) { st.time = c; st.timeUnit = Number(v.slice(5)); if (st.value === c) st.value = c === 1 ? 2 : 1; }
        else { if (st.time === c) st.time = 0; }
        if (st.value > s.nCols) st.value = 1;
        render();
      } else if (t.id === "imp-units") { st.scaling = Number(t.value); render(); }
    };
    dlg.oninput = (e) => {
      if (e.target.id === "imp-tau0") {
        st.tau0In = e.target.value;
        const { tau0 } = derived();
        dlg.querySelector("[data-act=load]").disabled = !(tau0 > 0);
      }
    };
    dlg.oncancel = (e) => { e.preventDefault(); finish(null); };

    render();
    dlg.showModal();
  });
}
