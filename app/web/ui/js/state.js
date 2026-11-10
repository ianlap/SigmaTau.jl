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
  band: true,              // intervals as a shaded band (else error bars)
  noiseId: true,           // noise-type strip and slope guides on the Stability plot
  fit: true,               // power-law model lines and the noise model table on Stability
  fitEdit: false,          // the noise model box shows the component editor
  mask: "",                // ITU-T wander mask drawn over MTIE / TDEV on Stability ("" = none)
  bias: "exact",           // bias correction of the totals and Thêo1, one of BIASES
  noise: "auto",           // noise type behind the intervals, one of NOISES
  tauMin: null,            // averaging-time limits in seconds, null = none
  tauMax: null,
  dataDrag: "zoom",        // Data plot drag: "zoom" | "range" (pick the analysis range)
  pins: [],                // pinned runs on Stability: { id, label, clocks: [{ name, slot, results }] }
  optionsOpen: false,      // options section when several clocks are loaded
  axis: "sigma_y",         // "sigma_y" | "sigma_x"
  tableDev: null,          // estimator in the per-τ table
  compareTau: null,        // τ of the "Compare at τ" table
  // Spectrum page (display settings, applied without a run)
  spec: { carrier: 10e6, carrierText: "10e6", phase: "L", avg: "log", guides: true },
  // Dynamic page: estimator, window = record / frac, map or 3D waterfall
  // Ensemble page: timescale settings and the last result ({ key, ts })
  ens: { algorithm: "at1", wmax: 1, pivot: 1, steps: "", result: null, error: null },
  dyn: { dev: "dadev", frac: 8, view: "map", cell: null, sliceTau: null },  // sliceTau: τ of the multi-clock comparison
  // Holdover page: training stretch = record / frac, line or parabola (drift), start from the last
  // reading (anchor), the time error limit in seconds (null = none) and the share of trials within it
  hold: { frac: 4, kalman: false, drift: false, anchor: true, limit: null, limitText: "", confidence: 0.95 },
  // Live page: where the readings come from and what they are
  live: { source: "sim", address: "192.168.1.20:5025", device: "/dev/ttyUSB0", query: "READ?", deviceQuery: "", reading: "hz",
          nominal: "", tau0: 0.1, devs: ["adev", "mdev"] },
  // Synthesis page: the record and the noise of each clock to be simulated (R, q1, q2, q3 are the
  // clock-model coefficients), and the predicted stability ({ key, doc }). The two clocks are the ones
  // the Ensemble page simulates on its own.
  synth: { n: 16384, tau0: 1, seed: 42, dev: "mdev", preview: null, error: null,
    clocks: [0, 1].map((j) => ({ name: `Sim ${j + 1}`, R: 0, q1: 1e-24 * 4 ** j, q2: 1e-30 * 4 ** j, q3: 0 })) },
  dataMode: "file",        // "file" | "synth" on the Data page
  dataView: "phase",       // "phase" | "drift": the Data page's plot
  driftModel: "linear",    // frequency drift model fitted on the drift view
  driftResid: false,       // drift view shows the residuals instead of the fit
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
    drift: "none",           // drift fitted and removed: "none" | "linear" | "quadratic"
    fillgaps: false,
    range: null,             // [start, end] in seconds from the first sample, or null = all
    remove_outliers: false,
    record: null,            // record + check from the core
    result: null,            // last analysis from the core (with the spectrum)
    dynamic: null,           // last dynamic-deviation map, with its settings key
    holdover: null,          // last holdover test, with its settings key
    driftFit: null,          // last frequency drift fit, with its settings key
    components: [],          // hand-built noise model: [{ noise, from, to }] in order; empty = automatic fit
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
  s.run.gen = (s.run.gen ?? 0) + 1;      // a run started before this is out of date
}
