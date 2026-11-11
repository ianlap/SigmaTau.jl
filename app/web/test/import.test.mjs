// What the import dialog guesses from a file's top: node --test app/web/test/import.test.mjs
import test from "node:test";
import assert from "node:assert/strict";

globalThis.window = {};
const { sniff, guessRoles, clockStats, timeStats, parseFactor, scaleFields } = await import("../ui/js/import.js");

const seq = (n, f) => Array.from({ length: n }, (_, i) => f(i));
const read = (name, lines) => { const s = sniff(name, lines.join("\n") + "\n"); return { s, g: guessRoles(s) }; };

test("a clock column that only drifts upward is not the time column", () => {
  const { s, g } = read("sats.csv", ["MJD,E08 phase (s),G08 phase (s)", ...seq(30, (i) => `${60904 + i / 2880},${1e-4 + i * 1e-9},${5e-4 + i * 2e-9}`)]);
  assert.equal(s.nCols, 3);
  assert.deepEqual([g.time, g.value, g.timeUnit], [1, 2, 86400]);
});

test("every column named for a reading is proposed as a clock, of the first one's kind", () => {
  const six = read("ens.csv", ["t (s),maser1 phase (s),maser2 phase (s),cs1 phase (s)", ...seq(30, (i) => `${i},${1e-9 * i},${2e-9 * i * i},${3e-9 * Math.sin(i)}`)]);
  assert.deepEqual([six.g.time, six.g.value, six.g.values], [1, 2, [2, 3, 4]]);
  const utc = read("labs.csv", ["MJD,UTC-UTC(NIST) (ns),UTC-UTC(PTB) (ns)", ...seq(30, (i) => `${50814 + 5 * i},${(i * 7) % 5},${(i * 3) % 4}`)]);
  assert.deepEqual([utc.g.values, utc.g.scaling], [[2, 3], 1e-9]);
  const mixed = read("m.csv", ["t,A phase (s),B freq,C phase (s)", ...seq(30, (i) => `${i},${1e-9 * i},${1e-12 * (i % 3)},${2e-9 * i}`)]);
  assert.deepEqual(mixed.g.values, [2, 4]);
  const flag = read("f.dat", ["# t  ratio  flag", ...seq(30, (i) => `${59646 + i / 86400}\t${(i % 3 - 1) * 1e-15}\t1`)]);
  assert.deepEqual(flag.g.values, [2]);
  assert.deepEqual(read("plain.csv", seq(20, (i) => `${i},${1e-9 * Math.cos(i)},${3e-9 * Math.sin(i)}`)).g.values, [2]);
});

test("MJD in a column headed t is days; seconds of day at 10 Hz stay seconds", () => {
  const mjd = read("r.dat", ["# t  ratio  flag", ...seq(30, (i) => `${(59638.801389 + i / 86400).toFixed(6)} ${-7e-15 + i * 1e-17} 1`)]);
  assert.deepEqual([mjd.g.time, mjd.g.timeUnit], [1, 86400]);
  const sod = read("s.csv", ["t,phase", ...seq(30, (i) => `${(43200 + i / 10).toFixed(1)},${1e-9 * i}`)]);
  assert.deepEqual([sod.g.time, sod.g.timeUnit], [1, 1]);
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
  assert.equal(parseFactor(" 2 * (3 + 4) / 7 - .5"), 1.5);
  assert.equal(parseFactor("500000"), 5e5);
  assert.ok(Math.abs(parseFactor("1/(2*pi*518.295836590863e12)") / 3.0709e-16 - 1) < 1e-4);
  for (const bad of ["", "-1", "0", "abc", "process.exit()", "1/0", "2**", "1/*2*/", "2**3", "(1", "1)", "1 2", "1e"]) assert.ok(Number.isNaN(parseFactor(bad)), bad);
});

test("reopening a clock with a typed scale keeps it in the factor field", () => {
  assert.deepEqual(scaleFields("phase", 1e-9), { scaling: 1e-9, factor: "" });
  assert.deepEqual(scaleFields("phase", 3.07e-16), { scaling: 1, factor: "3.07e-16" });
  assert.deepEqual(scaleFields("frequency", 1e-9, "5e5"), { scaling: 1e-9, factor: "5e5" });
  assert.equal(parseFactor(scaleFields("phase", 3.07e-16).factor), 3.07e-16);
});

test("a BIPM UTC - UTC(k) file is recognized by its header despite the '+/-' columns", () => {
  const rows = seq(30, (i) => `${50814 + 5 * i}   ${i % 4}.5   +/-   0.2   2.1   2.1   `);
  const s = sniff("utc-lab", ["  MJD  [UTC-UTC(LAB )]/ns  uA/ns  uB/ns  u/ns  \r", "\r", ...rows, "\x1a"].join("\n"));
  assert.equal(s.format, "bipm_utc");
  assert.deepEqual([s.nCols, s.nRows, s.names[1], s.kindHint, s.bad], [2, 30, "UTC-UTC(LAB) (ns)", "phase", null]);
  assert.deepEqual([guessRoles(s).time, guessRoles(s).value, guessRoles(s).timeUnit, guessRoles(s).scaling], [1, 2, 86400, 1e-9]);
});

test("a RINEX clock file lists its clocks and takes the clock's own epochs", () => {
  const rec = (kind, n, min, sec, x) => `${kind} ${n} 2025 08 17 00 ${String(min).padStart(2, "0")} ${String(sec).padStart(2, "0")}.000000  2   ${x}E-07  2.5E-11`;
  const lines = ["     3.04                 C                    M                   RINEX VERSION / TYPE", "                                                            END OF HEADER"];
  for (let k = 0; k < 8; k++) {
    const m = k >> 1, sec = (k & 1) * 30;
    lines.push(rec("AS", "G08", m, sec, 1 + k));
    if (!sec && m % 2 === 0) lines.push(rec("AR", "ALGO00CAN", m, 0, 2 + m));
  }
  const { s, g } = read("clk.clk", lines);
  assert.equal(s.format, "rinex_clk");
  assert.deepEqual(s.names, ["MJD", "ALGO00CAN", "G08"]);
  assert.equal(s.nRows, 8);
  assert.deepEqual(s.clocks, ["station · 2 epochs", "satellite · 8 epochs"]);
  assert.deepEqual([g.time, g.timeUnit], [1, 86400]);
  assert.equal(s.columns[0][1] - s.columns[0][0] > 0, true);
  assert.equal(s.columns[2][5], 6e-7);
  // Each clock on its own grid: the station every 2 min, the satellite every 30 s.
  assert.deepEqual(clockStats(s, 1, [2, 3], 86400).map((p) => [p.ts.tau0, p.n]), [[120, 2], [30, 8]]);
});

test("exact timestamps keep their own step on a short record", () => {
  assert.equal(timeStats(seq(10, (i) => 86400 * i)).tau0, 86400);       // not 86000
  assert.equal(timeStats([0.89, 4.03, 7.18]).tau0, 3.14);                  // ties as the reader rounds
  assert.equal(timeStats(seq(30, (i) => 1 + i / 86400).map((x) => Math.round(x * 1e6) / 1e6), 86400).tau0, 1);
});

test("a RINEX clock bias with a D exponent is read", () => {
  const { s } = read("old.clk", ["     2.00           C                   G                   RINEX VERSION / TYPE", "                                                            END OF HEADER",
    "AS G01  2025 08 17 23 59 30.000000  2   -.123000000000D-03  .250000000000D-10"]);
  assert.equal(s.columns[1][0], -1.23e-4);
});

test("a session saved with the on/off bias setting opens with its menu value", async () => {
  const { applySession } = await import("../ui/js/session.js");
  const s = {};
  applySession(s, { settings: { bias: false } });
  assert.equal(s.bias, "none");
});
