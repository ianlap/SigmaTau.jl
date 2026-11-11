// "Add estimators" dialog: search, grouped checkboxes, Done. Changes apply
// when the dialog closes.
import { GROUPS, ESTIMATORS } from "./estimators.js";
import { state, update, markStale } from "./state.js";

export function openPicker() {
  const dlg = document.getElementById("picker-dialog");
  const chosen = new Set(state.devs);

  const groupHtml = (g, q) => {
    const items = ESTIMATORS.filter((e) => e.group === g.id &&
      (!q || `${e.label} ${e.note} ${e.family} ${g.title}`.toLowerCase().includes(q)));
    if (!items.length) return "";
    return `<div class="picker-group"><h3>${g.title} · <i>${g.sub}</i></h3>${items.map((e) =>
      `<label class="pick" title="${e.family}"><input type="checkbox" value="${e.id}" ${chosen.has(e.id) ? "checked" : ""}><b>${e.label}</b><span>${e.note}</span></label>`).join("")}</div>`;
  };
  const body = (q = "") => GROUPS.map((g) => groupHtml(g, q.trim().toLowerCase())).join("") ||
    `<p class="muted">No estimator matches.</p>`;

  dlg.innerHTML = `
    <div class="picker-head">
      <div class="side-head"><h2>Add estimators</h2><span class="note">The choice applies when the dialog closes</span></div>
      <input type="text" id="picker-q" placeholder="Search, e.g. total, time, theo" autocomplete="off">
    </div>
    <div class="picker-body" id="picker-body">${body()}</div>
    <div class="picker-foot"><span class="muted small" id="picker-count"></span><button class="btn btn-primary" type="button" id="picker-done">Done</button></div>`;

  const count = () => {
    const slow = ESTIMATORS.filter((e) => e.slow && chosen.has(e.id)).map((e) => e.label);
    dlg.querySelector("#picker-count").textContent =
      `${chosen.size} selected${slow.length ? ` · ${new Intl.ListFormat("en-GB").format(slow)} take${slow.length > 1 ? "" : "s"} seconds on a long record` : ""}`;
  };
  count();

  const close = () => {
    dlg.close();
    const devs = ESTIMATORS.map((e) => e.id).filter((id) => chosen.has(id));
    // Keep the order the user already had, append new picks in catalog order.
    const order = [...state.devs.filter((d) => chosen.has(d)), ...devs.filter((d) => !state.devs.includes(d))];
    if (order.join() !== state.devs.join()) update((s) => { s.devs = order; markStale(s, "inputs"); });
  };

  dlg.querySelector("#picker-q").oninput = (e) => {
    dlg.querySelector("#picker-body").innerHTML = body(e.target.value);
  };
  dlg.onchange = (e) => {
    if (e.target.type !== "checkbox") return;
    e.target.checked ? chosen.add(e.target.value) : chosen.delete(e.target.value);
    count();
  };
  dlg.querySelector("#picker-done").onclick = close;
  dlg.oncancel = (e) => { e.preventDefault(); close(); };
  dlg.onclick = (e) => { if (e.target === dlg) close(); };
  dlg.showModal();
  dlg.querySelector("#picker-q").focus();
}
