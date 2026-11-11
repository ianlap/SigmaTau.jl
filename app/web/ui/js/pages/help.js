// Help page: the manual of the explorer, a topic at a time, at the address
// #help/topic/section. The text is help-content.js (see ../help.js).
import { TOPICS } from "../help-content.js";
import { find, address } from "../help.js";
import { escapeHtml as h } from "../format.js";

let root, topic;

export function mount(el) {
  root = el;
  topic = null;
  root.innerHTML = `<nav class="help-nav card" aria-label="Help topics"></nav><article class="help card"></article>`;
}

// Show the topic of the route, at its section.
export function route(parts) {
  const at = find(parts);
  if (at.topic !== topic) {
    topic = at.topic;
    const i = TOPICS.indexOf(topic);
    const turn = (t, text) => (t ? `<a href="${address(t)}">${text(t.nav)}</a>` : "<span></span>");
    root.querySelector(".help").innerHTML = `
      ${topic.sections.map((s) => `<section data-id="${h(s.id)}">${s.html.join("")}</section>`).join("")}
      <footer>${turn(TOPICS[i - 1], (t) => `← ${t}`)}
        <a href="${topic.web}" target="_blank" rel="noopener">This page in the documentation</a>
        ${turn(TOPICS[i + 1], (t) => `${t} →`)}</footer>`;
    root.querySelector(".help-nav").innerHTML = TOPICS.map((t) => `<a href="${address(t)}" class="${t === topic ? "on" : ""}">${t.nav}</a>${
      t === topic ? t.sections.slice(1).map((s) => `<a href="${h(address(t, s))}" class="sub l${s.level}">${s.title}</a>`).join("") : ""}`).join("");
  }
  const el = at.section.id && [...root.querySelectorAll("section")].find((s) => s.dataset.id === at.section.id);
  if (el) el.scrollIntoView(); else root.scrollTo(0, 0);
}

export function refresh() {}
export function hint() { return "The manual of the explorer. ? opens the help menu on any page"; }
