// The help: its topics against the pages, its addresses and its search.
//   node --test app/web/test/help.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";

const UI = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../ui");
const { TOPICS } = await import("../ui/js/help-content.js");
const { find, address, search, GREEK } = await import("../ui/js/help.js");

const text = TOPICS.flatMap((t) => t.sections.flatMap((s) => s.html)).join("\n");
const sources = ["index.html", ...readdirSync(path.join(UI, "js"), { recursive: true })
  .filter((f) => f.endsWith(".js") && f !== "help-content.js").map((f) => path.join("js", f))]
  .map((f) => readFileSync(path.join(UI, f), "utf8"));
const titles = (hits) => hits.map((h) => `${h.topic.nav} / ${h.section.title}`);

test("every page of the explorer has a topic, which F1 opens", () => {
  const pages = [...sources[0].matchAll(/data-page="(\w+)"/g)].map((m) => m[1]);
  assert.equal(pages.length, 9);
  for (const page of pages) assert.equal(find([page]).topic.id, page);
});

test("every address of the help, in its text and in the pages, names a section", () => {
  const used = [text, ...sources].flatMap((s) => [...s.matchAll(/href="(#help[^"]*)"/g)].map((m) => m[1]));
  for (const hash of ["#help", "#help/index/Keyboard-and-mouse", "#help/spectrum/Records-with-gaps"]) assert.ok(used.includes(hash), hash);
  for (const hash of new Set(used)) {
    const at = find(hash.split("/").slice(1));
    assert.equal(hash === "#help" ? "#help" : address(at.topic, at.section), hash);
  }
});

test("search: each word begins a word, accents aside, and a Greek letter answers to its name", () => {
  const grid = search("tau grid");
  assert.ok(titles(grid).includes("The side panel / Options"));
  assert.equal(grid[0].snippet[1], "τ");
  assert.ok(search("theo1").some((h) => h.snippet[1] === "Thêo1"));
  assert.ok(titles(search("export csv")).includes("Exports, sessions and reports / The report"));
  assert.deepEqual([search(""), search("zzz"), search("grid zzz")], [[], [], []]);
  assert.ok(search("(").length > 0);
  const greek = new Set(text.toLowerCase().match(/\p{Script=Greek}/gu));
  assert.deepEqual([...greek].sort(), Object.values(GREEK).sort());
});

test("search: the sections with the words in their title come first", () => {
  assert.deepEqual(titles(search("noise identification")).slice(0, 1), ["Stability / Noise identification"]);
  const hits = search("gaps");
  assert.ok(hits.length > 3 && hits[0].inTitle === 1 && hits.at(-1).inTitle === 0);
  assert.ok(hits.every((h, i) => !i || h.inTitle < hits[i - 1].inTitle || h.inTitle === hits[i - 1].inTitle && h.mentions <= hits[i - 1].mentions));
});
