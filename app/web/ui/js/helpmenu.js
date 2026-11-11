// The help menu under "?": a search through the manual, the sections on the
// page in view, and the ways into the help page.
import { TOPICS, SITE } from "./help-content.js";
import { search, address } from "./help.js";
import { state } from "./state.js";
import { escapeHtml as h } from "./format.js";

const btn = document.getElementById("help-btn");
const menu = document.getElementById("help-menu");
const input = document.getElementById("help-search");
const list = document.getElementById("help-items");

function items() {
  const query = input.value.trim();
  if (query) {
    const hits = search(query);
    if (!hits.length) return `<p class="muted">Nothing in the help matches.</p>`;
    return hits.map(({ topic, section, snippet: [lead, match, rest] }) => `
      <a href="${h(address(topic, section))}"><b>${section.title}</b>${section.id ? ` <span class="muted">${topic.nav}</span>` : ""}
        <small>${h(lead)}<mark>${h(match)}</mark>${h(rest)}</small></a>`).join("");
  }
  const here = TOPICS.find((t) => t.id === state.page);
  return `
    ${here ? `<h3>On this page</h3>${here.sections.filter((s) => s.level < 3).map((s) =>
      `<a href="${h(address(here, s))}">${s.id ? s.title : `About ${here.nav}`}</a>`).join("")}<hr>` : ""}
    <a href="#help">Contents of the help</a>
    <a href="#help/index/Keyboard-and-mouse">Keyboard and mouse</a>
    <a href="#help/estimators">Estimators and noise types</a>
    <hr>
    <a href="${SITE}" target="_blank" rel="noopener">Documentation on the web</a>`;
}

function show(open) {
  menu.hidden = !open;
  btn.setAttribute("aria-expanded", String(open));
  if (!open) return;
  input.value = "";
  list.innerHTML = items();
  input.focus();
}
export const openHelpMenu = () => show(true);

btn.addEventListener("click", () => show(menu.hidden));
input.addEventListener("input", () => { list.innerHTML = items(); });
menu.addEventListener("click", (e) => { if (e.target.closest("a")) show(false); });
menu.addEventListener("keydown", (e) => {
  if (e.key === "Escape") {
    e.stopPropagation();        // the key closes the menu and leaves a run alone
    show(false);
    btn.focus();
  } else if (e.key === "Enter" && e.target === input) {
    list.querySelector("a")?.click();
  } else if (e.key === "ArrowDown" || e.key === "ArrowUp") {
    e.preventDefault();
    const stops = [input, ...list.querySelectorAll("a")];
    stops[stops.indexOf(document.activeElement) + (e.key === "ArrowDown" ? 1 : -1)]?.focus();
  }
});
