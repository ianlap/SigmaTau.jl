// The frame around the page, plot first: a rail on the left with the
// panel's button, the results button and Run. The panel starts folded; a wide
// window can pin it beside the page (kept for the next visit), a narrower one
// opens it and the results card over the page, so the plot keeps the width.
import { state } from "./state.js";

const WIDE = matchMedia("(min-width: 1280px)");
const KEY = "sigmatau-panel";
const rail = document.getElementById("rail");
const btnPanel = document.getElementById("rail-panel");
const btnResults = document.getElementById("rail-results");
const btnRun = document.getElementById("rail-run");
let folded = true;
try { folded = localStorage.getItem(KEY) !== "pinned"; } catch { /* storage blocked: start folded */ }
let panelOpen = false, resultsOpen = false;

export function initLayout({ run }) {
  btnPanel.addEventListener("click", () => {
    if (WIDE.matches) {
      folded = !folded;
      try { localStorage.setItem(KEY, folded ? "folded" : "pinned"); } catch { /* not kept */ }
    } else {
      panelOpen = !panelOpen;
      resultsOpen = false;
    }
    apply();
  });
  btnResults.addEventListener("click", () => { resultsOpen = !resultsOpen; panelOpen = false; apply(); });
  btnRun.addEventListener("click", () => { closeOver(); run(); });
  WIDE.addEventListener("change", () => { panelOpen = resultsOpen = false; apply(); });
  // An open overlay closes with Esc or a click elsewhere (a dialog it opened excepted).
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && (panelOpen || resultsOpen) && state.run.state !== "running" && !document.querySelector("dialog[open]")) closeOver(); });
  document.addEventListener("pointerdown", (e) => {
    if (!(panelOpen || resultsOpen) || e.target.closest("#side, #rail, dialog, .page > .side-card")) return;
    closeOver();
  });
  apply();
}

// The focus, if it was inside, goes back to the rail button that opened it.
function closeOver() {
  if (!panelOpen && !resultsOpen) return;
  const back = document.activeElement?.closest("#side, .page > .side-card") && (panelOpen ? btnPanel : btnResults);
  panelOpen = resultsOpen = false;
  apply();
  back?.focus();
}

function apply() {
  const b = document.body, wide = WIDE.matches;
  const shown = wide ? !folded : panelOpen;
  b.classList.toggle("wide", wide);
  b.classList.toggle("panel-hidden", !shown);
  b.classList.toggle("panel-over", !wide && panelOpen);
  b.classList.toggle("results-over", !wide && resultsOpen);
  btnPanel.setAttribute("aria-expanded", String(shown));
  btnPanel.title = wide ? (folded ? "Pin the panel beside the page" : "Fold the panel away") : "Clocks and settings";
  btnResults.setAttribute("aria-expanded", String(resultsOpen));
  renderRail();
  window.dispatchEvent(new Event("resize"));   // the plots take the new width
}

// After each render: Run's state, and the results button only where the page has a results card to show.
export function renderRail() {
  const running = state.run.state === "running";
  btnRun.classList.toggle("running", running);
  btnRun.disabled = !running && !(state.clocks.length && state.devs.length);
  btnRun.title = running ? "Stop the analysis (F5 or Esc)" : "Run the analysis (F5)";
  btnRun.setAttribute("aria-label", running ? "Stop analysis" : "Run analysis");
  btnRun.querySelector(".stale-dot").hidden = running || !(state.run.stale && state.clocks.some((c) => c.record));
  btnResults.hidden = WIDE.matches || !document.querySelector(".page > .side-card:not([hidden]):not(.synth-clocks)");
  if (btnResults.hidden && resultsOpen) { resultsOpen = false; document.body.classList.remove("results-over"); }
  // A page with no panel, Run or results card on screen has no rail either.
  const empty = [btnPanel, btnResults, btnRun].every((b) => b.hidden || getComputedStyle(b).display === "none");
  if (rail.hidden !== empty) { rail.hidden = empty; requestAnimationFrame(() => window.dispatchEvent(new Event("resize"))); }
}
