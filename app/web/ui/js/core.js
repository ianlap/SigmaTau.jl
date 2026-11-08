// The only module that talks to the analysis core: serve.py over local HTTP
// in a browser, or the desktop shell's commands (app/desktop) inside its
// window. Both run sigmatau-cli and take and return the same shapes.
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

function recordRequest(clock) {
  return {
    source: clock.source.type === "synth"
      ? { type: "synth", n: clock.source.n, seed: clock.source.seed }
      : { type: "file", name: clock.source.name, content: clock.source.content },
    kind: clock.kind,
    tau0: clock.source.type === "synth" ? clock.source.tau0 : clock.tau0,
    time_col: clock.time_col,
    value_col: clock.value_col,
    reference_col: referenceColFor(clock),
    header: clock.header,
    scaling: clock.scaling,
    detrend: clock.detrend,
    fillgaps: clock.fillgaps,
    remove_outliers: clock.remove_outliers,
  };
}

// Read a record and its check (samples, span, gaps, outliers, drift),
// without running any estimator.
export function load(clock, signal) {
  return post({ ...recordRequest(clock), devs: [], series: 4000 }, signal);
}

// Run the selected estimators on one clock.
export function analyze(clock, signal) {
  return post({
    ...recordRequest(clock),
    devs: state.devs,
    grid: state.grid,
    confidence: state.ci ? state.confidence : null,
    bias_correction: state.bias,
    series: 0,
    spectrum: true,
  }, signal);
}

// Settings that decide a clock's dynamic map; a change means recompute.
export const dynamicKey = (clock) =>
  JSON.stringify([recordRequest(clock), state.dyn.dev, state.dyn.frac, state.grid]);

// Sliding-window deviation map for one clock (window = record / frac).
export function dynamic(clock, signal) {
  const n = clock.record?.n ?? 0;
  return post({
    ...recordRequest(clock),
    devs: [],
    grid: state.grid,
    series: 0,
    dynamic: state.dyn.dev,
    window: Math.max(16, Math.floor(n / state.dyn.frac)),
  }, signal);
}
