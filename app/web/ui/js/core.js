// The only module that talks to the analysis core. Today that is serve.py
// over local HTTP (it runs sigmatau-cli); a Tauri shell replaces the fetch
// in post() with `invoke` on a sidecar command and keeps the same shapes.
import { state, referenceColFor } from "./state.js";

export async function health() {
  try {
    const r = await fetch("api/health");
    return r.ok ? await r.json() : { ok: false, core: null };
  } catch {
    return { ok: false, core: null };
  }
}

async function post(body, signal) {
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
