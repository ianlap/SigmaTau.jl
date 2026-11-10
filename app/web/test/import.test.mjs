// What the import dialog guesses from a file's top: node --test app/web/test/import.test.mjs
import test from "node:test";
import assert from "node:assert/strict";

globalThis.window = {};
const { sniff, guessRoles, parseFactor } = await import("../ui/js/import.js");

const seq = (n, f) => Array.from({ length: n }, (_, i) => f(i));
const read = (name, lines) => { const s = sniff(name, lines.join("\n") + "\n"); return { s, g: guessRoles(s) }; };

test("a clock column that only drifts upward is not the time column", () => {
  const { s, g } = read("sats.csv", ["MJD,E08 phase (s),G08 phase (s)", ...seq(30, (i) => `${60904 + i / 2880},${1e-4 + i * 1e-9},${5e-4 + i * 2e-9}`)]);
  assert.equal(s.nCols, 3);
  assert.deepEqual([g.time, g.value, g.timeUnit], [1, 2, 86400]);
});

test("a blank-separated file named .csv is split on blanks", () => {
  const { s, g } = read("phase.csv", seq(20, (i) => `${i} ${35 + i / 10}`));
  assert.equal(s.nCols, 2);
  assert.equal(s.delimName, "whitespace-separated");
  assert.deepEqual([g.time, g.value], [1, 2]);
});

test("a one-column file gets its sample interval from a 'data interval' note", () => {
  const { s, g } = read("p.txt", ["# data interval 1.0 s", ...seq(20, (i) => `${1e-9 * i}`)]);
  assert.equal(s.tau0Hint, 1);
  assert.equal(g.time, 0);
});

test("units come from the column name", () => {
  const { g } = read("utc.csv", ["MJD,UTC-UTC(NIST) (ns)", ...seq(20, (i) => `${50814 + 5 * i},${(i * 7) % 5}`)]);
  assert.equal(g.scaling, 1e-9);
});

test("readings far below a nanosecond are frequency, and a comment line can name the columns", () => {
  const { s, g } = read("ratio.dat", ["# t  ratio  flag", ...seq(20, (i) => `${59646 + i / 86400}\t${(i % 3 - 1) * 1e-15}\t1`)]);
  assert.deepEqual(s.names, ["t", "ratio", "flag"]);
  assert.equal(s.kindHint, "frequency");
  assert.deepEqual([g.time, g.value], [1, 2]);
});

test("a typed scale is a number or a short expression, nothing else", () => {
  assert.equal(parseFactor("1e-9"), 1e-9);
  assert.ok(Math.abs(parseFactor("1/(2*pi*518.295836590863e12)") / 3.0709e-16 - 1) < 1e-4);
  for (const bad of ["", "-1", "0", "abc", "process.exit()", "1/0", "2**"]) assert.ok(Number.isNaN(parseFactor(bad)), bad);
});
