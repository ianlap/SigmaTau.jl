// Application state: one plain object plus a change signal. Pages read
// `state` directly and call `update()` to change it; listeners re-render.

const listeners = new Set();

export const state = {
  page: "stability",
  core: { ok: false, name: null },
  clocks: [],              // see newClock()
  devs: ["adev", "mdev", "hdev", "mhdev"],
  grid: "octave",
  ci: true,
  confidence: 0.683,
  band: false,
  bias: true,
  axis: "sigma_y",         // "sigma_y" | "sigma_x"
  layout: "overlay",       // "overlay" | "side"
  tableClock: null,
  tableDev: null,
  run: { state: "idle", text: "", lastMs: null, autoRun: true },
};

let nextId = 1;
export function newClock(fields) {
  return {
    id: nextId++,
    name: `Clock ${String.fromCharCode(64 + nextId - 1)}`,
    visible: true,
    source: { type: "synth", n: 16384, seed: 42 },
    kind: "phase",
    tau0: null,           // null = from the file
    time_col: 1,
    value_col: 2,
    header: 0,
    scaling: 1,
    detrend: "none",
    result: null,         // last response from the core
    error: null,
    ...fields,
  };
}

export function subscribe(fn) {
  listeners.add(fn);
  return () => listeners.delete(fn);
}

// `what` tells listeners what changed: "analysis" means inputs to the
// core changed (so results are stale), "view" means only presentation.
export function update(mutator, what = "view") {
  mutator(state);
  for (const fn of listeners) fn(what);
}

export const activeClocks = () => state.clocks.filter((c) => c.visible);
