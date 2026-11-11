// Screenshots for the explorer manual (docs/src/explorer/), written to
// docs/src/assets/explorer/. Every record shown is simulated here, from fixed
// seeds.
//
//   SIGMATAU_CLI=cli/build/bin/sigmatau node docs/screenshots.mjs
//
// Needs Playwright with Chromium, as app/web/test/explorer.e2e.mjs does, and
// python3 with Pillow, which reduces each image to 256 colors.
import { spawn, execFileSync } from "node:child_process";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";
import path from "node:path";

const { chromium } = createRequire(import.meta.url)("playwright");
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const OUT = path.join(ROOT, "docs/src/assets/explorer");
const PORT = 9300 + Math.floor(Math.random() * 600);
const BASE = `http://127.0.0.1:${PORT}/`;

// ── Simulated records ──
// Normal deviates from a seed (mulberry32, Box–Muller).
function normals(seed) {
  let a = seed >>> 0;
  const uniform = () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  return () => Math.sqrt(-2 * Math.log(1 - uniform())) * Math.cos(2 * Math.PI * uniform());
}

// Phase in seconds at 1 s: white PM `pm` (s), white FM `fm` and random-walk FM
// `rw` (per sample), with the white FM multiplied by `step` over the last 40 %
// of the record.
function phase(n, seed, { pm = 0, fm = 1e-12, rw = 0, step = 1 } = {}) {
  const g = normals(seed), x = [];
  let t = 0, walk = 0;
  for (let i = 0; i < n; i++) {
    walk += rw * g();
    t += fm * (i >= 0.6 * n ? step : 1) * g() + walk;
    x.push(t + pm * g());
  }
  return x;
}

const file = (name, text) => ({ name, mimeType: "text/plain", buffer: Buffer.from(text) });

const table = (x) => x.map((v, i) => `${i} ${v.toExponential(9)}`).join("\n") + "\n";

// One clock: three noise types and two outliers.
function clockFile() {
  const x = phase(32768, 1, { pm: 4e-12, fm: 2e-12, rw: 1e-14 });
  x[9000] += 3e-9; x[21000] -= 3e-9;
  return file("clock.dat", "# time_s phase_s\n" + table(x));
}

// One clock whose white FM triples at 60 % of the record.
const stepFile = () => file("step.dat", table(phase(16384, 2, { step: 3 })));

// Three clocks against one reference.
function threeFile() {
  const cols = [phase(20000, 3, { fm: 1e-12, rw: 2e-15 }), phase(20000, 4, { fm: 2e-12, rw: 1e-15 }), phase(20000, 5, { fm: 5e-13, rw: 4e-15 })];
  const rows = cols[0].map((_, i) => [i, ...cols.map((c) => c[i].toExponential(9))].join(","));
  return file("three_clocks.csv", `time_s,clock_A,clock_B,clock_C\n${rows.join("\n")}\n`);
}

// ── The explorer ──
const server = spawn("python3", [path.join(ROOT, "app/web/serve.py"), "--port", String(PORT)], { stdio: "ignore" });
const browser = await chromium.launch();
const made = [];
try {
  for (let i = 0; i < 100; i++) {
    try { if ((await (await fetch(BASE + "api/health")).json()).ok) break; } catch { /* not up yet */ }
    await new Promise((r) => setTimeout(r, 100));
  }
  const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1.5, colorScheme: "light" });
  await context.addInitScript(() => localStorage.clear());   // no import settings remembered from scene to scene
  const page = await context.newPage();
  const SLOW = { timeout: 300_000 };   // an uncompiled core takes seconds per call
  const settled = async () => {
    await page.waitForTimeout(200);
    await page.waitForFunction(() => !/Reading|Running|Computing|Starting/.test(document.getElementById("status-text").textContent), null, SLOW);
    await page.waitForTimeout(900);   // plots drawn, fonts in
  };
  const plotted = (sel) => page.waitForFunction((s) => document.querySelector(s)?.data?.length > 0, sel, SLOW);
  const start = async (name) => { await page.goto("about:blank"); await page.goto(`${BASE}#${name}`); await settled(); };
  const go = async (name) => { await page.click(`a[data-page=${name}]`); await settled(); };
  const open = async (f) => {
    const [chooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("button[data-act=open]")]);
    await chooser.setFiles(f);
    await page.waitForSelector("#import-dialog[open]");
  };
  const load = async () => { await page.click("#import-dialog [data-act=load]"); await settled(); };
  const run = async () => { await page.click(".run-btn"); await settled(); };
  const shot = async (name, target = page) => {
    await page.mouse.move(0, 0);   // no hover label, no plot toolbar
    await page.waitForTimeout(300);
    made.push(path.join(OUT, `${name}.png`));
    await target.screenshot({ path: made.at(-1) });
  };

  // The window at start.
  await start("stability");
  await shot("window");

  // Synthesis: white PM added to the first clock.
  await go("synthesis");
  await page.fill(".synth-clock[data-i='0'] [data-q=R]", "1e-23");
  await page.locator(":focus").blur();
  await settled();
  await shot("synthesis");

  // Ensemble, on its own four simulated clocks.
  await start("ensemble");
  await plotted("[data-ref=adev]");
  await settled();
  await shot("ensemble");

  // Live: the simulator, half a minute at ten readings a second.
  await go("live");
  await page.click("[data-ctl=go]");
  await page.waitForTimeout(30_000);
  await shot("live");
  await page.click("[data-ctl=go]");

  // Data, Stability and Spectrum on one clock.
  await start("data");
  await open(clockFile());
  await load();
  await shot("data", page.locator("#page"));
  await page.click("[data-act=remove-outliers]");
  await run();
  await go("stability");
  await shot("stability", page.locator("#page"));
  await go("spectrum");
  await shot("spectrum");

  // Holdover on the same clock, with a limit of 1 ns.
  await go("holdover");
  await page.fill("[data-ctl=limit]", "1e-9");
  await page.locator(":focus").blur();
  await page.waitForFunction(() => document.querySelector("#page [data-ref=empty]")?.hidden
    && document.querySelector("#page [data-ref=plot]")?.data?.length > 0
    && /Within \d/.test(document.querySelector("#page [data-ref=side]")?.textContent ?? ""), null, SLOW);
  await settled();
  await shot("holdover");

  // Steering on the same clock.
  await go("steering");
  await page.waitForFunction(() => document.querySelector("#page [data-ref=empty]")?.hidden
    && /Time constant/.test(document.querySelector("#page [data-ref=side]")?.textContent ?? ""), null, SLOW);
  await settled();
  await shot("steering");

  // Dynamic: a clock whose noise changes during the record.
  await start("data");
  await open(stepFile());
  await load();
  await run();
  await go("dynamic");
  await plotted(".dyn-single [data-ref=map]");
  await settled();
  const map = await page.locator(".dyn-single [data-ref=map] .nsewdrag").boundingBox();
  await page.mouse.click(map.x + 0.3 * map.width, map.y + 0.8 * map.height);
  await settled();
  await shot("dynamic");

  // Three clocks from one file: the import dialog and the comparison.
  await start("data");
  await open(threeFile());
  for (const col of [3, 4]) await page.selectOption(`#import-dialog select[data-col="${col}"]`, "value");
  await page.fill("#imp-against", "reference");
  await page.locator(":focus").blur();
  await shot("import", page.locator("#import-dialog"));
  await load();
  for (const dev of ["hdev", "mhdev"]) await page.click(`[data-remove=${dev}]`);
  await run();
  await go("stability");
  await shot("stability-clocks");
} finally {
  await browser.close();
  server.kill();
}

// 256 colors by median cut: a good third of the size, and the colors of the
// page surfaces stay as they are.
execFileSync("python3", ["-c", `
import sys
from PIL import Image
for f in sys.argv[1:]:
    im = Image.open(f).convert("RGB")
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(f, optimize=True)
`, ...made]);
console.log(`${made.length} screenshots in ${path.relative(ROOT, OUT)}`);
