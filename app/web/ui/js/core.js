// The only module that talks to the analysis core: serve.py over local HTTP
// in a browser, or the desktop shell's commands (app/desktop) inside its
// window. Both run the sigmatau binary (cli/) with the arguments built here
// and hand back its JSON.
import { state, referenceColFor } from "./state.js";

const tauri = window.__TAURI__?.core;
let nextRun = 1;

// One analysis through the desktop shell; aborting stops the core process.
async function invokeAnalyze(body, signal) {
  const id = nextRun++;
  const abortError = () => new DOMException("The analysis was cancelled.", "AbortError");
  if (signal?.aborted) throw abortError();
  const cancelled = new Promise((_, reject) => signal?.addEventListener("abort", () => {
    tauri.invoke("cancel", { id }).catch(() => {});
    reject(abortError());
  }, { once: true }));
  const out = await Promise.race([tauri.invoke("analyze", { id, req: body }), cancelled]);
  if (out.error) throw new Error(out.error);
  return out;
}

export async function health() {
  if (tauri) {
    try { return await tauri.invoke("health"); } catch { return { ok: false, core: null }; }
  }
  try {
    const r = await fetch("api/health");
    return r.ok ? await r.json() : { ok: false, core: null };
  } catch {
    return { ok: false, core: null };
  }
}

async function post(body, signal) {
  if (tauri) {
    try {
      return await invokeAnalyze(body, signal);
    } catch (err) {
      throw err instanceof Error ? err : new Error(String(err));
    }
  }
  const r = await fetch("api/analyze", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
    signal,
  });
  let out;
  try {
    out = await r.json();
  } catch {
    throw new Error(`the server answered ${r.status}`);
  }
  if (!r.ok || out.error) throw new Error(out.error || `the server answered ${r.status}`);
  return out;
}

// Command-line flags that describe one clock's record. A file travels
// separately (the shell saves it and passes its path as FILE).
function recordArgs(clock) {
  const synth = clock.source.type === "synth";
  const a = synth ? ["--synth", clock.source.n, "--seed", clock.source.seed] : [];
  const opt = (flag, v) => { if (v !== null && v !== undefined && v !== "") a.push(flag, v); };
  opt("--tau0", synth ? clock.source.tau0 : clock.tau0);
  // An instrument file (TimeLab, Stable32 .phd/.frd, counter log) is read by
  // its layout, which fixes the record type, units and columns.
  if (clock.format) {
    a.push("--input-format", clock.format);
    // A table is read by column name, as the record type the user picked.
    if (clock.format === "table") {
      opt("--column", clock.column_name);
      a.push(clock.kind === "frequency" ? "--freq" : "--phase");
    }
  } else recordColumns(clock, a, opt);
  if (clock.drift && clock.drift !== "none") a.push("--remove-drift", clock.drift);
  if (clock.fillgaps) a.push("--fillgaps");
  if (clock.remove_outliers) a.push("--remove-outliers");
  a.push("--format", "json");
  return a;
}

function recordColumns(clock, a, opt) {
  if (clock.kind === "frequency") a.push("--freq");
  opt("--time-col", clock.time_col);
  opt("--value-col", clock.value_col);
  const ref = referenceColFor(clock);
  if (ref) opt("--reference-col", ref);
  opt("--header", clock.header);
  opt("--scaling", clock.scaling);
}

// One run of the sigmatau binary: a command, its flags, and the record file.
function request(clock, command, extra) {
  const file = clock.source.type === "synth" ? null : { name: clock.source.name, content: clock.source.content };
  return { args: [command, ...recordArgs(clock), ...extra].map(String), file };
}

// The core answers in its versioned JSON (schema_version 1, see cli/README.md).
// The pages read columns, so it is reshaped once, here; `doc` keeps the
// original for the JSON export.
function view(doc) {
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  const r = doc.record ?? {};
  const record = {
    n: r.n, tau0: r.tau0, kind: r.quantity, source: doc.source, series: r.series,
    check: {
      samples: r.n, span_s: r.span_s, gaps: r.gaps, outliers: r.outliers,
      outlier_t: r.outlier_t ?? [], outlier_x: r.outlier_x ?? [], drift_per_day: r.drift_per_day,
      removed_offset: r.removed_offset ?? null, removed_drift: r.removed_drift_per_day ?? null,
    },
  };
  const results = (doc.results ?? []).map((res) => {
    const col = (k) => res.points.map((p) => p[k]);
    // Columns an estimator does not report (no CI for MTIE and TIE rms, or
    // --no-ci) are all null; the pages expect them empty.
    const has = (k) => res.points.some((p) => p[k] !== null && p[k] !== undefined);
    const ci = has("ci_lower");
    return {
      dev: res.deviation, label: res.label, axis: res.axis, tau: col("tau"), sigma: col("dev"),
      ci_lo: ci ? col("ci_lower") : [], ci_hi: ci ? col("ci_upper") : [],
      edf: has("edf") ? col("edf") : [], noise: has("noise_type") ? col("noise_type") : [], neff: col("neff"),
    };
  });
  const d = doc.dynamic;
  const dynamic = d && { dev: d.deviation, window: d.window, step: d.step, t: d.t, tau: d.tau, map: d.dev };
  return { record, results, spectrum: doc.spectrum, dynamic, drift: doc.drift, elapsed_ms: doc.elapsed_ms, doc };
}

const run = async (body, signal) => view(await post(body, signal));

// Read a record and its check (samples, span, gaps, outliers, drift),
// without running any estimator.
export function load(clock, signal) {
  return run(request(clock, "check", ["--series", 4000]), signal);
}

// Run the selected estimators on one clock, with its spectrum.
export function analyze(clock, signal) {
  return run(request(clock, "analyze", [
    "--devs", state.devs.join(","),
    "--grid", state.grid,
    ...(state.ci ? ["--confidence", state.confidence] : ["--no-ci"]),
    ...(state.bias ? [] : ["--no-bias-correction"]),
    "--spectrum",
  ]), signal);
}

// Fit a frequency drift model (linear, quadratic, log, diffusion) to one clock.
export function drift(clock, signal) {
  return run(request(clock, "drift", ["--model", state.driftModel, "--series", 2000]), signal);
}

// Settings that decide a clock's drift fit; a change means refit.
export const driftKey = (clock) => JSON.stringify([request(clock, "drift", []), state.driftModel]);

// Settings that decide a clock's dynamic map; a change means recompute.
export const dynamicKey = (clock) =>
  JSON.stringify([request(clock, "dynamic", []), state.dyn.dev, state.dyn.frac, state.grid]);

// Sliding-window deviation map for one clock (window = record / frac).
export function dynamic(clock, signal) {
  const n = clock.record?.n ?? 0;
  return run(request(clock, "dynamic", [
    "--grid", state.grid,
    "--dev", state.dyn.dev,
    "--window", Math.max(16, Math.floor(n / state.dyn.frac)),
  ]), signal);
}

// ── Live measurement ──
// `sigmatau live` runs until stopped and writes one JSON line per update;
// the shell keeps the lines and the page polls for the ones it has not seen.
export async function liveStart(args) {
  if (tauri) return tauri.invoke("live_start", { args });
  const r = await fetch("api/live/start", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ args }) });
  const out = await r.json();
  if (!r.ok) throw new Error(out.error || `the server answered ${r.status}`);
  return out.id;
}

export async function livePoll(id, since) {
  const out = tauri ? await tauri.invoke("live_poll", { id, since })
    : await (await fetch(`api/live/poll?id=${id}&since=${since}`)).json();
  const lines = (out.lines ?? []).map((l) => JSON.parse(l));
  const failed = lines.find((l) => l.error);
  return { lines: lines.filter((l) => !l.error), next: out.next, running: out.running, error: failed?.error ?? out.error ?? null };
}

export function liveStop() {
  return tauri ? tauri.invoke("live_stop")
    : fetch("api/live/stop", { method: "POST", headers: { "Content-Type": "application/json" }, body: "{}" });
}
