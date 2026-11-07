// Application state: one plain object plus a change signal. Code reads
// `state` directly and changes it through update(); listeners re-render.

const listeners = new Set();

export const state = {
  page: "stability",
  core: { ok: false, name: null },
  clocks: [],              // see newClock()
  activeClock: null,       // id of the clock the Data page shows
  reference: null,         // id of the reference clock, or null (as measured)
  devs: ["adev", "mdev", "hdev", "mhdev"],
  grid: "octave",
  ci: true,
  confidence: 0.683,
  band: false,
  bias: true,
  optionsOpen: false,      // options section when several clocks are loaded
  axis: "sigma_y",         // "sigma_y" | "sigma_x"
  tableDev: null,          // estimator in the per-τ table
  compareTau: null,        // τ of the "Compare at τ" table
  dataMode: "file",        // "file" | "synth" on the Data page
  run: { state: "idle", stale: null, note: "" },  // stale: null | "record" | "inputs"
};

let nextId = 1;
export function newClock(fields = {}) {
  const used = new Set(state.clocks.map((c) => c.slot));
  let slot = 0;
  while (used.has(slot)) slot++;
  return {
    id: nextId++,
    slot,                    // colour index: --s1 … --s8, then repeats
    name: "Clock",
    visible: true,
    source: { type: "synth", n: 16384, seed: 42, tau0: 1 },
    sourceInfo: "",          // one line under the file name
    kind: "phase",
    tau0: null,              // null = from the file
    time_col: 1,
    value_col: 2,
    group: null,             // clocks read from one file share a group
    against: "",             // what the file says the clocks were measured against
    header: 0,
    scaling: 1,
    detrend: "none",
    fillgaps: false,
    remove_outliers: false,
    hasTimeCol: false,
    record: null,            // record + check from the core
    result: null,            // last analysis from the core
    error: null,
    ...fields,
  };
}

export function subscribe(fn) {
  listeners.add(fn);
  return () => listeners.delete(fn);
}

export function update(mutator) {
  mutator(state);
  for (const fn of listeners) fn();
}

export const referenceClock = () => state.clocks.find((c) => c.id === state.reference) ?? null;
export const isReference = (c) => c.id === state.reference;

// Column to subtract from this clock's value column: the reference clock's,
// when both come from the same file. 0 = none.
export function referenceColFor(c) {
  const r = referenceClock();
  if (!r || r === c || c.source.type !== "file" || r.group === null || r.group !== c.group) return 0;
  return r.value_col;
}

// Clocks that get analysed and plotted: not hidden, not the reference.
export const visibleClocks = () => state.clocks.filter((c) => c.visible && !isReference(c));
export const multi = () => state.clocks.length > 1;
export const activeClock = () =>
  state.clocks.find((c) => c.id === state.activeClock) ?? state.clocks[0] ?? null;

// Mark results out of date. "record" (a new or changed record) outranks
// "inputs" (estimators or options).
export function markStale(s, why) {
  if (s.run.stale !== "record") s.run.stale = why;
}
