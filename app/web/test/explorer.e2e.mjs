// End-to-end smoke test of the explorer in a browser, against a real analysis
// core: load a fixture on the Data page, run the analysis and check the
// Stability plot against the command line's own JSON, then visit Spectrum,
// Dynamic, Live (simulated source) and Help, come back to Stability, and run
// a record with gaps, with no page errors.
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

  await t.test("Live: a simulated source streams and stops", async () => {
    await page.click("a[data-page=live]");
    await page.click("[data-ctl=go]");
    await plotted("#page [data-ref=sigma]");
    await page.click("[data-ctl=go]");
    await page.waitForFunction(() => /Stopped/.test(document.querySelector("#page .live-state")?.textContent ?? ""), null, SLOW);
  });

  await t.test("Help opens", async () => {
    await page.click("#help-btn");
    await page.waitForFunction(() => document.querySelector("#page")?.className.includes("page-help"));
    assert.ok((await page.locator("#page").innerText()).length > 200);
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
    await page.locator("#import-dialog .gap-card", { hasText: "1 zero reading" }).waitFor();
    await page.fill("#imp-tau0", "1");
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /z\.frd.*51[23] samples/ }).waitFor();   // the new record, not the last one
    assert.ok(await page.locator("input[data-field=fillgaps]").isChecked());
  });

  await t.test("Data: remembered .frd settings still show a file's zeros", async () => {
    const y = Array.from({ length: 256 }, (_, i) => (i % 100 === 50 ? "0" : String(1e-12 * (2 + Math.cos(i / 5)))));
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "z2.frd", mimeType: "text/plain", buffer: Buffer.from(y.join("\n")) });
    await page.locator("#import-dialog .gap-card", { hasText: "3 zero readings" }).waitFor();
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /z2\.frd.*25[67] samples/ }).waitFor();
  });

  await t.test("Data: a one-column frequency file with NaN readings is filled", async () => {
    await page.goto(base + "#data");
    const y = Array.from({ length: 512 }, (_, i) => (i === 100 ? "NaN" : String(1e-12 * Math.sin(i / 3))));
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles({ name: "y.txt", mimeType: "text/plain", buffer: Buffer.from(y.join("\n")) });
    await page.fill("#imp-tau0", "1");
    await page.getByRole("button", { name: "Load record" }).click();
    await settled();
    await page.locator("#record-chip-text", { hasText: /y\.txt.*51[23] samples/ }).waitFor();   // the new record, not the last one
    assert.ok(await page.locator("input[data-field=fillgaps]").isChecked());
  });

  assert.deepEqual(errors, []);
});
