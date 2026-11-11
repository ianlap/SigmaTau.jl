// The help: its topics and sections, and the search through them. The text is
// the manual of the explorer. docs/explorer_help.jl writes help-content.js
// from docs/src/explorer, a topic for each page and a section for each heading.
import { TOPICS } from "./help-content.js";

// The topic and section of the route "help/topic/section", as its parts after
// "help". An unknown topic is the overview, an unknown section the top.
export function find([id, ...rest] = []) {
  const topic = TOPICS.find((t) => t.id === id) ?? TOPICS[0];
  const section = topic.sections.find((s) => s.id === rest.join("/")) ?? topic.sections[0];
  return { topic, section };
}

export const address = (topic, section) => `#help/${topic.id}${section?.id ? `/${section.id}` : ""}`;

// ── Search ──
const ENTITY = { amp: "&", lt: "<", gt: ">", quot: '"' };
const plain = (html) => html.replace(/<\/(p|li|td|th|h\d|pre)>/g, " ").replace(/<[^>]+>/g, "")
  .replace(/&(amp|lt|gt|quot);/g, (_, e) => ENTITY[e]).replace(/\s+/g, " ").trim();
// Text as it is compared: lower case and without accents, each character still in its place.
const fold = (s) => s.normalize("NFD").replace(/\p{M}/gu, "").toLowerCase();
// A Greek letter of the manual answers to its name: "tau grid" finds the τ grid.
// A subscript follows its letter in the text: "sigma_y" and "σ_y" find σy.
export const GREEK = { alpha: "α", chi: "χ", nu: "ν", omega: "ω", pi: "π", sigma: "σ", tau: "τ" };

let index = null;      // every section with its text, made at the first search

// The sections in which every word of the query begins a word of the section,
// its title included, or of the name of its page, and one word at least of
// the section itself. Those with more of the words in their title come first,
// then those whose text has the words more often. Each comes with the text
// around its first match, as [before, match, after].
export function search(query) {
  const words = fold(query).split(/\s+/).filter(Boolean).map((w) => {
    const forms = [...new Set([w, w.replaceAll("_", ""), w.split("_").map((p) => GREEK[p] ?? p).join("")])]
      .filter(Boolean).map((f) => f.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"));
    return new RegExp(`(^|[^\\p{L}\\p{N}])(${forms.join("|")})`, "gu");
  });
  if (!words.length) return [];
  index ??= TOPICS.flatMap((topic) => topic.sections.map((section) => {
    const text = plain(section.html.slice(1).join(" ")).normalize("NFC");
    return { topic, section, text, title: fold(plain(section.title)), page: fold(topic.nav), body: fold(text) };
  }));
  const hits = [];
  for (const { topic, section, text, title, page, body } of index) {
    const found = words.map((w) => [...body.matchAll(w)]);
    const inTitle = words.filter((w) => title.search(w) >= 0).length;
    const mentions = found.reduce((n, f) => n + f.length, 0);
    if (!(inTitle || mentions)) continue;
    if (!words.every((w, i) => found[i].length || title.search(w) >= 0 || page.search(w) >= 0)) continue;
    const m = found.find((f) => f.length)?.[0];
    const from = m ? m.index + m[1].length : 0, to = from + (m ? m[2].length : 0);
    const lead = text.slice(Math.max(0, from - 48), from).replace(/^\S*\s/, from > 48 ? "… " : "$&");
    const snippet = [lead, text.slice(from, to), text.slice(to, to + 96) + (text.length > to + 96 ? "…" : "")];
    hits.push({ topic, section, snippet, inTitle, mentions });
  }
  return hits.sort((a, b) => b.inTitle - a.inTitle || b.mentions - a.mentions).slice(0, 30);
}
