// The only module that talks to the analysis core. Today that is serve.py
// over local HTTP (it runs sigmatau-cli); a Tauri shell would replace the two
// fetch calls with `invoke` on a sidecar command and keep the same shapes.
import { state } from "./state.js";

export async function health() {
  try {
    const r = await fetch("api/health");
    if (!r.ok) return { ok: false, core: null };
    return await r.json();
  } catch {
    return { ok: false, core: null };
  }
}

export function requestFor(clock) {
  return {
    source: clock.source,
    kind: clock.kind,
    tau0: clock.tau0,
    time_col: clock.time_col,
    value_col: clock.value_col,
    header: clock.header,
    scaling: clock.scaling,
    detrend: clock.detrend,
    devs: state.devs,
    grid: state.grid,
    confidence: state.ci ? state.confidence : null,
    bias_correction: state.bias,
    series: 4000,
  };
}

export async function analyze(clock, signal) {
  const r = await fetch("api/analyze", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(requestFor(clock)),
    signal,
  });
  let body;
  try {
    body = await r.json();
  } catch {
    throw new Error(`server answered ${r.status}`);
  }
  if (!r.ok || body.error) throw new Error(body.error || `server answered ${r.status}`);
  return body;
}
