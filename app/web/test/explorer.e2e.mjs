// End-to-end smoke test of the explorer in a browser, against a real analysis
// core: load a fixture on the Data page, run the analysis and check the
// Stability plot and the Holdover and Steering curves against the command
// line's own JSON, then visit Spectrum, Dynamic, Ensemble (a Kalman scale with
// flicker states, then the hat against the command line), Live (simulated
// source) and Help, come back to Stability for a pinned run and a telecom
// mask, and run a record and a file of three clocks with gaps, with no page
// errors.
//
//   SIGMATAU_CLI=cli/build/bin/sigmatau node --test app/web/test/explorer.e2e.mjs
//
// Needs Playwright with Chromium (npm install playwright; npx playwright install chromium).
import test from "node:test";
import assert from "node:assert/strict";
import { spawn, execFileSync } from "node:child_process";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";
import path from "node:path";

const { chromium } = createRequire(import.meta.url)("playwright");
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../../..");
const FIXTURE = path.join(ROOT, "test/fixtures/validation/stable32gen.DAT");
const CLI = process.env.SIGMATAU_CLI;
const PORT = 9300 + Math.floor(Math.random() * 600);
const SLOW = { timeout: 300_000 };   // an uncompiled core takes seconds per call

test("explorer pages, against the command line", { skip: !CLI && "set SIGMATAU_CLI", timeout: 900_000 }, async (t) => {
  const server = spawn("python3", [path.join(ROOT, "app/web/serve.py"), "--port", String(PORT)], { stdio: "ignore" });
  t.after(() => server.kill());
  const base = `http://127.0.0.1:${PORT}/`;
  for (let i = 0; i < 100; i++) {
    try { if ((await (await fetch(base + "api/health")).json()).ok) break; } catch { /* not up yet */ }
    await new Promise((r) => setTimeout(r, 100));
  }
  const browser = await chromium.launch();
  t.after(() => browser.close());
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  // Web fonts, and the vendored Plotly that falls back to the CDN when ui/vendor is empty, may 404.
  page.on("console", (m) => m.type() === "error" && !/fonts\.g|vendor\/plotly/.test(m.location().url) && errors.push(m.text()));
  const settled = () => page.waitForFunction(() => !/Reading|Running|Computing|Starting/.test(document.getElementById("status-text").textContent), null, SLOW);
  const plotted = (sel) => page.waitForFunction((s) => document.querySelector(s)?.data?.length > 0, sel, SLOW);

  await t.test("Data: a Stable32 file through the import dialog", async () => {
    await page.goto(base + "#data");
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles(FIXTURE);
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    assert.match(await page.locator("#record-chip-text").innerText(), /stable32gen\.DAT.*8,?\s?192 samples/);
  });

  await t.test("Stability: the plotted σ equal the command line's", async () => {
    await page.click("a[data-page=stability]");
    await page.click(".run-btn");
    await plotted("#page .js-plotly-plot");
    await settled();
    const shown = await page.evaluate(() => Object.fromEntries(document.querySelector("#page .js-plotly-plot").data
      .filter((d) => /^[A-Z]+$/.test(d.name ?? "")).map((d) => [d.name.toLowerCase(), { tau: d.x, dev: d.y }])));
    const doc = JSON.parse(execFileSync(CLI, ["analyze", FIXTURE, "--devs", "adev,mdev,hdev,mhdev", "--format", "json"]));
    assert.deepEqual(Object.keys(shown).sort(), ["adev", "hdev", "mdev", "mhdev"]);
    for (const r of doc.results) {
      assert.deepEqual(shown[r.deviation].tau, r.points.map((p) => p.tau), r.deviation);
      assert.deepEqual(shown[r.deviation].dev, r.points.map((p) => p.dev), r.deviation);
    }
  });

  await t.test("Spectrum and Dynamic draw their plots", async () => {
    await page.click("a[data-page=spectrum]");
    await plotted("#page .js-plotly-plot");
    await page.click("a[data-page=dynamic]");
    await plotted("#page [data-ref=map]");
  });

  await t.test("Holdover: the plotted RMS equals the command line's", async () => {
    await page.click("a[data-page=holdover]");
    await plotted("#page [data-ref=plot]");
    const rms = await page.evaluate(() => document.querySelector("#page [data-ref=plot]").data.find((d) => d.name === "RMS").y);
    const doc = JSON.parse(execFileSync(CLI, ["holdover", FIXTURE, "--format", "json"]));
    assert.deepEqual(rms, doc.holdover.rms);
    // A stretch of twice the holdover time: equal to √6·TDEV at τ₀.
    await page.click("[data-ctl=frac][data-v='0']");
    await page.waitForFunction(() => document.querySelector("#page .hold-card .cell-kv")?.innerText.includes("2τ"));
    const at0 = await page.evaluate(() => ["RMS", "√6·TDEV"].map((n) => document.querySelector("#page [data-ref=plot]").data.find((d) => d.name === n).y[0]));
    const grown = JSON.parse(execFileSync(CLI, ["holdover", FIXTURE, "--grow", "2", "--format", "json"])).holdover;
    assert.deepEqual(at0, [grown.rms[0], grown.reference[0]]);
    assert.ok(Math.abs(at0[0] / at0[1] - 1) < 1e-9);
    // A Kalman filter as the predictor: its claim equals the command line's.
    await page.click('[data-ctl=kalman][data-v="1"]');
    // The curve is dotted past the longest τ of the fit its model comes from.
    const claimed = await (await page.waitForFunction(() => {
      const d = document.querySelector("#page [data-ref=plot]").data, fit = d.find((t) => t.name === "Claimed by the filter")?.y;
      return fit && [fit, d.find((t) => t.name === "Claimed, past the fit").y.slice(1)];
    }, null, SLOW)).jsonValue();
    const kalman = JSON.parse(execFileSync(CLI, ["holdover", FIXTURE, "--kalman", "--format", "json"])).holdover;
    assert.deepEqual(claimed.flat(), kalman.claimed);
    assert.equal(claimed[0].length, kalman.tau.filter((t) => t <= kalman.fit_tau_max).length);
  });

  await t.test("Steering: the plotted σ equal the command line's", async () => {
    await page.click("a[data-page=steering]");
    await plotted("#page.page-steering [data-ref=plot]");
    const curve = (name) => page.evaluate((n) => document.querySelector("#page [data-ref=plot]").data.find((d) => d.name === n)?.y, name);
    const args = ["steer", FIXTURE, "--format", "json"];   // the page's stretch, a quarter of the record, is the default
    const doc = JSON.parse(execFileSync(CLI, [...args, "--time-constant", "64"])).steer;
    assert.deepEqual(await curve("Steered"), doc.steered);
    assert.deepEqual(await curve("Predicted (clock model)"), doc.predicted.adev);
    // Gains from a cost: the loop card shows the command line's phase gain.
    await page.click('[data-ctl=from][data-v="cost"]');
    await page.waitForFunction(() => /Rings/.test(document.querySelector("#page .steer-card")?.innerText ?? ""), null, SLOW);
    const cost = JSON.parse(execFileSync(CLI, [...args, "--cost", "1:0:1e8"])).steer;
    assert.deepEqual(await curve("Steered"), cost.steered);
    assert.ok((await page.locator("#page .steer-card").innerText()).includes(cost.gains.phase.toExponential(3)));
    // The steered phase over the record, with the correction on its own axis.
    await page.click('[data-ctl=view][data-v="record"]');
    await page.waitForFunction(() => document.querySelector("#page [data-ref=plot]").data.some((d) => d.yaxis === "y2"));
  });

  await t.test("Ensemble: a Kalman scale with flicker states", async () => {
    await page.click("a[data-page=ensemble]");
    await page.selectOption("#side [data-ctl=alg]", "greenhall");
    await page.selectOption("#side [data-ctl=flicker]", "on");
    await page.waitForFunction(() => /Greenhall.*flicker states/.test(document.querySelector("#page [data-ref=sub]")?.textContent ?? ""), null, SLOW);
    assert.match(await page.locator("#page .ens-note").innerText(), /[Ff]licker/);
    assert.match(await page.locator("#page [data-ref=table]").innerText(), /White FM/);   // the fitted models
    await plotted("#page [data-ref=adev]");
  });

  await t.test("Ensemble: the hat's σ of each clock equal the command line's", async () => {
    await page.click("[data-ctl=view][data-v=hat]");
    await page.waitForFunction(() => document.querySelector("#page [data-ref=adev]").data?.some((d) => d.name === "Pairs, as measured"), null, SLOW);
    assert.equal(await page.locator("#page [data-ref=title]").innerText(), "ADEV of each clock");
    const shown = await page.evaluate(() => document.querySelector("#page [data-ref=adev]").data.filter((d) => /^clock \d$/.test(d.name)).map((d) => d.y));
    const { hat } = JSON.parse(execFileSync(CLI, ["hat", "--synth", "16384", "--seed", "42", "--format", "json"]));   // the page's own four clocks
    assert.deepEqual(shown, hat.results[0].clocks.map((c) => c.points.map((p) => p.dev)));
    await page.click("[data-ctl=view][data-v=scale]");
  });

  await t.test("Live: a simulated source streams and stops", async () => {
    await page.click("a[data-page=live]");
    await page.click("[data-ctl=go]");
    await plotted("#page [data-ref=sigma]");
    await page.click("[data-ctl=go]");
    await page.waitForFunction(() => /Stopped/.test(document.querySelector("#page .live-state")?.textContent ?? ""), null, SLOW);
  });

  await t.test("Help: the menu lists the page in view, a search opens a section, F1 the page's own help", async () => {
    await page.keyboard.press("?");
    await page.waitForSelector("#help-menu:not([hidden])");
    assert.match(await page.locator("#help-items").textContent(), /On this page\s*About Live[\s\S]*While it runs/);
    await page.fill("#help-search", "pinned runs");
    await page.keyboard.press("Enter");
    await page.waitForFunction(() => location.hash === "#help/stability/Pinned-runs" && document.querySelector("#page").className.includes("page-help"));
    assert.ok(await page.locator("#help-menu").isHidden());
    const at = await page.evaluate(() => ({
      heading: document.querySelector('#page section[data-id="Pinned-runs"] h2').getBoundingClientRect().top,
      button: document.getElementById("help-btn").getBoundingClientRect().right,
    }));
    const { width, height } = page.viewportSize();
    assert.ok(at.heading >= 0 && at.heading < height, `the section is in view (its heading at ${at.heading} px)`);
    assert.ok(at.button <= width, `the ? button is inside the window (its right edge at ${at.button} px)`);
    await page.click("a[data-page=spectrum]");
    await page.waitForSelector("#page.page-spectrum");
    await page.keyboard.press("F1");
    await page.waitForFunction(() => location.hash === "#help/spectrum");
    assert.deepEqual(await page.locator("#page article h1, #page .help-nav a.on").allInnerTexts(), ["Spectrum", "Spectrum"]);
    // An address that names no section of the topic is put right.
    await page.evaluate(() => { location.hash = "#help/spectrum/nowhere"; });
    await page.waitForFunction(() => location.hash === "#help/spectrum");
  });

  await t.test("Stability, revisited: one click pins one run and opens the model editor", async () => {
    await page.click("a[data-page=stability]");
    await page.click("[data-act=pin]");
    assert.equal(await page.locator("#page [data-unpin]").count(), 1);
    await page.click("#page [data-fit=toggle]");
    assert.equal(await page.locator("#page .model-edit").count(), 1);
    await page.click("#page [data-fit=toggle]");
  });

  await t.test("Stability: a telecom mask over MTIE, with the command line's verdict", async () => {
    await page.click("[data-act=pick]");
    await page.check("#picker-dialog input[value=mtie]");
    await page.click("#picker-done");
    await page.click("[data-axis=sigma_x]");
    await page.selectOption("[data-ref=mask]", "G8272_PRTC_A");
    await page.click(".run-btn");
    const card = page.locator("#page .model-box", { hasText: "PRTC-A" });
    await card.waitFor(SLOW);
    const mask = JSON.parse(execFileSync(CLI, ["analyze", FIXTURE, "--devs", "mtie", "--mask", "G8272_PRTC_A", "--format", "json"])).masks[0];
    const text = await card.innerText();
    assert.ok(text.includes(mask.source) && text.includes(mask.checks[0].pass ? "passes" : "fails"), text);
    assert.ok(await page.evaluate(() => document.querySelector("#page .js-plotly-plot").data.some((d) => d.name === "G.8272 PRTC-A MTIE mask")));
    await page.click("[data-axis=sigma_y]");
  });

  await t.test("Stability: a phase-noise table is drawn as the command line converts it", async () => {
    const table = path.join(mkdtempSync(path.join(tmpdir(), "sigmatau-")), "ocxo.csv");
    writeFileSync(table, "Offset (Hz),L(f) (dBc/Hz)\n1,-100\n10,-130\n100,-150\n1000,-160\n1e4,-165\n1e5,-165\n");
    await page.setInputFiles("#page [data-ref=pn-file]", table);
    await page.waitForFunction(() => document.querySelectorAll("#page [data-unpin]").length === 2, null, SLOW);
    assert.match(await page.locator("#page .chip").last().getAttribute("title"), /from 1 Hz to 100 kHz, carrier 10 MHz[^]*hold for 100 kHz/);
    const shown = await page.evaluate(() => document.querySelector("#page .js-plotly-plot").data.find((d) => d.name === "ℒ(f) 2 · MDEV"));
    const doc = JSON.parse(execFileSync(CLI, ["phasenoise", table, "--header", "1", "--carrier", "10e6", "--devs", "mdev", "--format", "json"]));
    assert.deepEqual([shown.x, shown.y], [doc.results[0].tau, doc.results[0].dev]);
    assert.equal(await page.evaluate(() => document.querySelector("#page .js-plotly-plot").layout.xaxis.dtick), 2, "nine decades of τ, a label every other one");
    await page.click("#page .chip:last-child [data-unpin]");
    // A table of one decade gives one τ, drawn as a point.
    writeFileSync(table, "1,-100\n10,-130\n");
    await page.setInputFiles("#page [data-ref=pn-file]", table);
    await page.waitForFunction(() => document.querySelectorAll("#page [data-unpin]").length === 2, null, SLOW);
    const point = await page.evaluate(() => document.querySelector("#page .js-plotly-plot").data.find((d) => d.name === "ℒ(f) 3 · MDEV"));
    assert.deepEqual([point.x.length, point.mode], [1, "markers"]);
    await page.click("#page .chip:last-child [data-unpin]");
  });

  await t.test("A record with gaps runs without its spectrum", async () => {
    const gapped = path.join(mkdtempSync(path.join(tmpdir(), "sigmatau-")), "gapped.csv");
    const rows = Array.from({ length: 1200 }, (_, i) => i).filter((i) => i < 500 || i >= 520);
    writeFileSync(gapped, rows.map((i) => `${i},${(1e-9 * Math.sin(i * 0.7) + 1e-12 * i * i).toExponential(6)}`).join("\n") + "\n");
    await page.click("a[data-page=data]");
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles(gapped);
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.click("a[data-page=stability]");
    await page.click(".run-btn");
    await settled();
    assert.match(await page.locator("#status-text").innerText(), /Up to date/);
    await page.click("a[data-page=spectrum]");
    await page.waitForFunction(() => /record has gaps/.test(document.querySelector("#page [data-ref=empty]")?.textContent), null, SLOW);
    // The core refuses the map and names a command-line flag; the page names the box that sets it.
    const n = errors.length;
    await page.click("a[data-page=dynamic]");
    await page.locator("#page.page-dynamic", { hasText: "without gaps; turn on Fill gaps by reflection on the Data page" }).waitFor(SLOW);
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
  });

  await t.test("A failed run: its error stays in the status bar until the box it names is ticked", async () => {
    await page.click("a[data-page=stability]");
    await page.click("[data-act=pick]");
    await page.check("#picker-dialog input[value=totdev]");   // cannot skip gaps
    await page.click("#picker-done");
    const n = errors.length;
    await page.click(".run-btn");
    const failed = page.locator("#status-text", { hasText: "turn on Fill gaps by reflection on the Data page" });
    await failed.waitFor(SLOW);
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
    await page.click("a[data-page=data]");
    assert.equal(await failed.count(), 1);
    await page.click("input[data-field=fillgaps]");
    await settled();
    assert.match(await page.locator("#status-hint").innerText(), /^Run analysis \(F5\)/);
  });

  await t.test("Data: a zero in a Stable32 frequency file is a gap, filled", async () => {
    await page.goto(base + "#data");
    const y = Array.from({ length: 512 }, (_, i) => (i === 100 ? "0" : String(1e-12 * (2 + Math.sin(i / 3)))));
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "z.frd", mimeType: "text/plain", buffer: Buffer.from(y.join("\n")) });
    await page.locator("#import-dialog .gap-card", { hasText: "1 zero reading" }).waitFor(SLOW);
    await page.fill("#imp-tau0", "1");
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /z\.frd.*51[23] samples/ }).waitFor(SLOW);   // the new record, not the last one
    assert.ok(await page.locator("input[data-field=fillgaps]").isChecked());
  });

  await t.test("Data: remembered .frd settings still show a file's zeros", async () => {
    const y = Array.from({ length: 256 }, (_, i) => (i % 100 === 50 ? "0" : String(1e-12 * (2 + Math.cos(i / 5)))));
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "z2.frd", mimeType: "text/plain", buffer: Buffer.from(y.join("\n")) });
    await page.locator("#import-dialog .gap-card", { hasText: "3 zero readings" }).waitFor(SLOW);
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /z2\.frd.*25[67] samples/ }).waitFor(SLOW);
  });

  await t.test("Data: a one-column frequency file with NaN readings is filled", async () => {
    await page.goto(base + "#data");
    const y = Array.from({ length: 512 }, (_, i) => (i === 100 ? "NaN" : String(1e-12 * Math.sin(i / 3))));
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "y.txt", mimeType: "text/plain", buffer: Buffer.from(y.join("\n")) });
    await page.fill("#imp-tau0", "1");
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /y\.txt.*51[23] samples/ }).waitFor(SLOW);   // the new record, not the last one
    assert.ok(await page.locator("input[data-field=fillgaps]").isChecked());
  });

  await t.test("Data: Fill gaps turned off and on again; the next file starts without it", async () => {
    const n = errors.length;
    await page.uncheck("input[data-field=fillgaps]");
    await page.locator(".check-card", { hasText: "could not read" }).waitFor(SLOW);
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));   // the core refuses the unfilled record
    await page.click("input[data-field=fillgaps]");   // (click: the card re-renders as the record reloads)
    await page.locator("#record-chip-text", { hasText: /y\.txt.*51[23] samples/ }).waitFor(SLOW);
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    const rows = Array.from({ length: 600 }, (_, i) => i).filter((i) => i < 300 || i >= 310);
    // (read with the settings remembered from gapped.csv, without the dialog)
    await chooser.setFiles({ name: "p.csv", mimeType: "text/csv", buffer: Buffer.from(rows.map((i) => `${i},${(1e-9 * Math.sin(i * 0.7)).toExponential(6)}`).join("\n")) });
    await settled();
    await page.locator("#record-chip-text", { hasText: /p\.csv/ }).waitFor(SLOW);
    assert.ok(!(await page.locator("input[data-field=fillgaps]").isChecked()));
  });

  await t.test("Ensemble: a file with gaps forms a Kalman scale and an averaging scale, and its clocks are separated", async () => {
    const rows = Array.from({ length: 1200 }, (_, i) => i).filter((i) => i < 500 || i >= 520);
    const clock = (i, k) => (1e-9 * Math.sin(i * (0.7 + 0.2 * k)) + 1e-12 * k * i).toExponential(6);
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "three.csv", mimeType: "text/csv", buffer: Buffer.from(rows.map((i) => [i, clock(i, 1), clock(i, 2), clock(i, 3)].join(",")).join("\n")) });
    for (const col of [3, 4]) await page.selectOption(`#import-dialog select[data-col="${col}"]`, "value");
    await page.click("#import-dialog [data-act=load]");
    await settled();
    const scale = (alg) => page.waitForFunction((a) => new RegExp(`^${a}.* 1200 samples`).test(document.querySelector("#page [data-ref=sub]")?.textContent ?? ""), alg, SLOW);
    await page.click("a[data-page=ensemble]");   // Greenhall, as the Ensemble step above left it
    await scale("Greenhall");
    await page.selectOption("#side [data-ctl=alg]", "at1");
    await scale("AT1");
    // The hat and the covariance separate the three clocks and the reference; the KLTS interval says why it cannot.
    const separated = (method) => page.waitForFunction((m) => document.querySelector("#page [data-ref=sub]")?.textContent.startsWith(m), method, SLOW);
    await page.click("[data-ctl=view][data-v=hat]");
    await separated("Hat · 4 clocks");
    await page.selectOption("#side [data-ctl=method]", "groslambert");
    await separated("Groslambert covariance · 4 clocks");
    const n = errors.length;
    await page.uncheck("#side [data-ctl=reference]");
    await page.selectOption("#side [data-ctl=method]", "klts");
    await page.locator("#page [data-ref=empty]", { hasText: "The KLTS interval needs a file without missing samples" }).waitFor(SLOW);
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
  });

  await t.test("Data at 960 px: a long file name and the plot's header stay inside their cards", async () => {
    await page.setViewportSize({ width: 960, height: 800 });
    await page.click("a[data-page=data]");
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    const row = (i) => [i, ...[1, 2, 3].map((k) => (1e-9 * Math.sin(0.7 * k * i)).toExponential(6))].join(",");
    // three clocks, read with the settings remembered from three.csv
    await chooser.setFiles({ name: "three_clocks_with_a_long_name.csv", mimeType: "text/csv", buffer: Buffer.from(Array.from({ length: 600 }, (_, i) => row(i)).join("\n")) });
    await settled();
    await page.locator("#side .clock-row .clock-text").last().click();   // a clock named after the file
    assert.match(await page.locator("[data-ref=plot-title]").innerText(), /three_clocks_with_a_long_name col 4$/);
    const outside = await page.evaluate(() => [...document.querySelectorAll("#page .card")].flatMap((card) => {
      const edge = card.getBoundingClientRect().right + 1;
      return [...card.querySelectorAll(".file-title *, .card-row *")].filter((e) => e.getBoundingClientRect().right > edge).map((e) => e.className || e.tagName);
    }));
    assert.deepEqual(outside, []);
    assert.ok((await page.locator("[data-ref=plot-hint]").boundingBox()).width >= 100, "the hint keeps 100 px");
  });

  await t.test("Data: a file of named clocks is proposed whole; Holdover and Steering name the clock they cannot compute", async () => {
    // Clock a has no readings in its first 800 samples, so no stretch to learn a noise model from.
    const clock = (i, k) => (1e-9 * Math.sin(i * (0.7 + 0.2 * k)) + 1e-12 * k * i).toExponential(6);
    const rows = Array.from({ length: 1200 }, (_, i) => [i, i < 800 ? "NaN" : clock(i, 1), clock(i, 2), clock(i, 3)].join(","));
    await page.goto(base + "#data");
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "late.csv", mimeType: "text/csv", buffer: Buffer.from(`t (s),a phase (s),b phase (s),c phase (s)\n${rows.join("\n")}\n`) });
    await page.getByRole("button", { name: "Load 3 clocks" }).click();   // all three, none marked by hand
    await settled();
    const n = errors.length;
    await page.click("a[data-page=holdover]");
    await page.click('[data-ctl=kalman][data-v="1"]');
    const note = page.locator("#page [data-ref=note]", { hasText: "a phase (s)" });
    await note.waitFor(SLOW);
    assert.match(await note.innerText(), /noise model has 0 usable τ points/);
    await page.waitForFunction(() => { const n = document.querySelector("#page [data-ref=plot]")?.data?.map((d) => d.name) ?? []; return ["b", "c"].every((k) => n.some((x) => x.startsWith(`${k} phase (s)`))); }, null, SLOW);
    assert.ok(!(await page.evaluate(() => document.querySelector("#page [data-ref=plot]").data.some((d) => d.name.startsWith("a phase (s)")))));
    await page.click("a[data-page=steering]");
    await page.locator("#page [data-ref=note]", { hasText: "a phase (s)" }).waitFor(SLOW);
    await plotted("#page [data-ref=plot]");
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
  });

  await t.test("Data: with no drift fit to show, the buttons above the plot stay clickable", async () => {
    await page.goto(base + "#data");
    const n = errors.length;
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "tiny.csv", mimeType: "text/csv", buffer: Buffer.from("0,1e-9\n1,3e-9\n2,2e-9\n") });
    await page.locator("#import-dialog [data-act=load]").click({ timeout: 5000 }).catch(() => {});   // (the dialog, unless settings were remembered)
    await settled();
    const reachable = () => page.evaluate(() => [...document.querySelectorAll("[data-view]")].map((b) => { const r = b.getBoundingClientRect(); return document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2)?.closest("[data-view]") === b; }));
    await page.click("[data-view=drift]");
    await page.locator("#page [data-ref=empty]", { hasText: "No fit for this record" }).waitFor(SLOW);
    for (const view of ["drift", "hist", "acf", "phase"]) {   // whatever the view shows, the view buttons stay reachable
      await page.click(`[data-view=${view}]`);
      await settled();
      assert.deepEqual(await reachable(), [true, true, true, true], view);
    }
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));   // the core refuses a record and a fit of three samples
  });

  await t.test("Live: the readings analyzed are the file that Save readings writes", async () => {
    await page.click("a[data-page=live]");
    const n = errors.length;
    await page.click("[data-ctl=go]");
    await plotted("#page [data-ref=sigma]");
    await page.waitForTimeout(3000);   // readings enough to analyze (16)
    await page.click("[data-ctl=go]");
    await page.locator("[data-act=analyze]:not([disabled])").waitFor(SLOW);
    const [download] = await Promise.all([page.waitForEvent("download"), page.click("[data-act=save]")]);
    const file = readFileSync(await download.path());
    const saved = createHash("sha256").update(file).digest("hex");
    const readings = file.toString().split("\n").filter((l) => l && !l.startsWith("#")).length;
    await page.click("[data-act=analyze]");
    await settled();
    const analyzed = await page.evaluate(async (url) => {
      const { state } = await import(url);
      const bytes = new TextEncoder().encode(state.clocks.at(-1).source.content);
      return [...new Uint8Array(await crypto.subtle.digest("SHA-256", bytes))].map((b) => b.toString(16).padStart(2, "0")).join("");
    }, base + "js/state.js");
    assert.equal(analyzed, saved);
    // Analyze readings runs every loaded clock: the 422s are the records the earlier steps left (three samples, TOTDEV on gaps), not the Live record.
    assert.equal(await page.evaluate(async (url) => (await import(url)).state.clocks.at(-1).record?.n, base + "js/state.js"), readings);
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
  });

  await t.test("A small window keeps the Data plot usable in every view and the hat's legend off the curves", async () => {
    await page.setViewportSize({ width: 960, height: 640 });
    await page.goto(base + "#data");
    const n = errors.length;
    let seed = 7;
    const noise = () => ((seed = (seed * 16807) % 2147483647) / 2147483647 - 0.5);
    const rows = Array.from({ length: 600 }, (_, i) => [i, ...[1, 2, 3].map(() => 1e-9 * noise())].join(",")).join("\n");
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "three_clocks_in_a_long_named_file.csv", mimeType: "text/csv", buffer: Buffer.from(rows + "\n") });
    await page.locator("#import-dialog [data-act=load]").click({ timeout: 5000 }).catch(() => {});
    await settled();
    for (const view of ["phase", "drift", "hist", "acf"]) {
      await page.click(`[data-view=${view}]`);
      await settled();
      const h = await page.evaluate(() => document.querySelector("#page [data-ref=plot]").getBoundingClientRect().height);
      assert.ok(h >= 250, `${view}: the plot is ${h} px high`);
    }
    await page.click("a[data-page=ensemble]");
    await page.click("[data-ctl=view][data-v=hat]");
    await page.selectOption("#side [data-ctl=method]", "hat");   // (an earlier step left KLTS, which a file with gaps refuses)
    await page.waitForFunction(() => document.querySelector("#page [data-ref=sub]")?.textContent.startsWith("Hat · "), null, SLOW);
    await plotted("#page [data-ref=adev]");
    const [legend, area] = await page.evaluate(() => ["#page [data-ref=adev] .legend", "#page [data-ref=adev] .nsewdrag"].map((q) => document.querySelector(q).getBoundingClientRect().toJSON()));
    assert.ok(legend.bottom <= area.top + 1, "the legend sits above the plot area");
    await page.setViewportSize({ width: 1280, height: 800 });
    assert.ok(errors.splice(n).every((e) => / 422 /.test(e)));
  });

  assert.deepEqual(errors, []);
});
