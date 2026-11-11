// End-to-end smoke test of the explorer in a browser, against a real analysis
// core: load a fixture on the Data page, run the analysis and check the
// Stability plot and the Holdover and Steering curves against the command
// line's own JSON, then visit Spectrum, Dynamic, Ensemble (a Kalman scale with flicker states),
// Live (simulated source) and Help, come back to Stability, and run a record
// with gaps, with no page errors.
//
//   SIGMATAU_CLI=cli/build/bin/sigmatau node --test app/web/test/explorer.e2e.mjs
//
// Needs Playwright with Chromium (npm install playwright; npx playwright install chromium).
import test from "node:test";
import assert from "node:assert/strict";
import { spawn, execFileSync } from "node:child_process";
import { mkdtempSync, writeFileSync } from "node:fs";
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
    const claimed = (await page.waitForFunction(() => document.querySelector("#page [data-ref=plot]").data.find((d) => d.name === "Claimed by the filter")?.y, null, SLOW)).jsonValue();
    assert.deepEqual(await claimed, JSON.parse(execFileSync(CLI, ["holdover", FIXTURE, "--kalman", "--format", "json"])).holdover.claimed);
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
    await plotted("#page [data-ref=adev]");
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

  assert.deepEqual(errors, []);
});
