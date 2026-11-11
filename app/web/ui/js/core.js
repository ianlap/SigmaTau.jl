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
  const abortError = () => new DOMException("The analysis was canceled.", "AbortError");
  if (signal?.aborted) throw abortError();
  const canceled = new Promise((_, reject) => signal?.addEventListener("abort", () => {
    tauri.invoke("cancel", { id }).catch(() => {});
    reject(abortError());
  }, { once: true }));
  const out = await Promise.race([tauri.invoke("analyze", { id, req: body }), canceled]);
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

// The flags that make a simulated record (or ensemble) of `src`: a
// "--noise-model R:q1:q2:q3" per clock when it was made on the Synthesis page.
export function synthArgs(src) {
  const a = ["--synth", src.n, "--seed", src.seed];
  if (src.models) for (const m of src.models) a.push("--noise-model", m.join(":"));
  else if (src.clocks > 1) a.push("--clocks", src.clocks);
  return a;
}

// Command-line flags that describe one clock's record. A file travels
// separately (the shell saves it and passes its path as FILE).
function recordArgs(clock) {
  const synth = clock.source.type === "synth";
  const a = synth ? synthArgs(clock.source) : [];
  // One clock of a simulated ensemble: the core simulates them all and hands back this one.
  if (synth && clock.source.clocks > 1) a.push("--clock", clock.source.index);
  const opt = (flag, v) => { if (v !== null && v !== undefined && v !== "") a.push(flag, v); };
  // (a RINEX clock file has a grid per clock: the reader takes it from the epochs)
  if (clock.format !== "rinex_clk") opt("--tau0", synth ? clock.source.tau0 : clock.tau0);
  // An instrument file (TimeLab, Stable32 .phd/.frd, counter log) is read by
  // its layout, which fixes the record type, units and columns.
  if (clock.format) {
    a.push("--input-format", clock.format);
    // A table is read by column name, as the record type the user picked;
    // a RINEX clock file by clock name.
    if (clock.format === "table") {
      opt("--column", clock.column_name);
      a.push(clock.kind === "frequency" ? "--freq" : "--phase");
    } else if (clock.format === "rinex_clk") opt("--column", clock.column_name);
  } else recordColumns(clock, a, opt);
  if (clock.drift && clock.drift !== "none") a.push("--remove-drift", clock.drift);
  if (clock.keep_zeros) a.push("--keep-zeros");
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

// The core answers in its versioned JSON (schema_version 1, see docs/src/reference/cli.md).
// The pages read columns, so it is reshaped once, here; `doc` keeps the
// original for the JSON export.
function view(doc) {
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  const r = doc.record ?? {};
  const record = {
    n: r.n, tau0: r.tau0, kind: r.quantity, source: doc.source, series: r.series,
    check: {
      samples: r.n, span_s: r.span_s, start: r.start || "", gaps: r.gaps, outliers: r.outliers,
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
  return { record, results, spectrum: doc.spectrum, powerlaw: doc.powerlaw, masks: doc.masks ?? [], dynamic, drift: doc.drift, elapsed_ms: doc.elapsed_ms, doc };
}

// For a record with gaps the core names the command line's remedy. A clock's
// request carries the preparation of the Data page, where the remedy is a box.
const GAP_REMEDY = /fill them with (--fillgaps|fillgaps first)$/;
const run = (body, signal) => post(body, signal).then(view, (err) => {
  throw GAP_REMEDY.test(err.message) ? new Error(err.message.replace(GAP_REMEDY, "turn on Fill gaps by reflection on the Data page")) : err;
});

// Read a record and its check (samples, span, gaps, outliers, drift),
// without running any estimator.
export function load(clock, signal) {
  return run(request(clock, "check", ["--series", 4000]), signal);
}

// --range (seconds from the first sample) for runs on part of a record. The
// record check reads the whole record, so the Data page can show the range.
const rangeArgs = (clock) => clock.range ? ["--range", `${clock.range[0]}:${clock.range[1]}`] : [];
const tauArgs = () => [
  ...(state.tauMin ? ["--tau-min", state.tauMin] : []),
  ...(state.tauMax ? ["--tau-max", state.tauMax] : []),
];

// Run the selected estimators on one clock, with its spectrum.
const analyzeRequest = (clock) => request(clock, "analyze", [
  ...rangeArgs(clock), ...tauArgs(),
  "--devs", state.devs.join(","),
  "--grid", state.grid,
  ...(state.ci ? ["--confidence", state.confidence] : ["--no-ci"]),
  "--bias", state.bias,
  "--noise", state.noise,
  "--spectrum", "--fit", "--mask", "all",
  ...clock.components.flatMap((k) => ["--component", `${k.noise}:${k.from ?? ""}:${k.to ?? ""}`]),
]);

export const analyze = (clock, signal) => run(analyzeRequest(clock), signal);

// The sigmatau command line that reproduces a clock's analysis (CSV out),
// for reports.
export function analyzeCommand(clock) {
  const { args, file } = analyzeRequest(clock);
  const valued = new Set(["--format", "--component", "--mask"]);   // JSON-only options and their values
  const kept = args.filter((a, i) => a !== "--spectrum" && a !== "--fit" && !valued.has(a) && !valued.has(args[i - 1]));
  const q = (a) => (/^[\w.,:+\/=-]+$/.test(a) ? a : `'${a.replace(/'/g, "'\\''")}'`);
  return ["sigmatau", kept[0], ...(file ? [file.name] : []), ...kept.slice(1)].map(q).join(" ");
}

// The estimator applies to the linear model; the others are least-squares fits.
export const driftMethod = () => (state.driftModel === "linear" ? state.driftMethod : "frequency_fit");

// Fit a frequency drift model (linear, quadratic, log, diffusion) to one clock.
export function drift(clock, signal) {
  return run(request(clock, "drift", ["--model", state.driftModel, "--method", driftMethod(), "--series", 2000]), signal);
}

// Settings that decide a clock's drift fit; a change means refit.
export const driftKey = (clock) => JSON.stringify([request(clock, "drift", []), state.driftModel, driftMethod()]);

// Settings that decide a clock's dynamic map; a change means recompute.
export const dynamicKey = (clock) =>
  JSON.stringify([request(clock, "dynamic", rangeArgs(clock)), state.dyn.dev, state.dyn.frac, state.grid]);

// Samples of the record, or of its analysis range.
const samples = (clock) =>
  clock.range ? Math.floor((clock.range[1] - clock.range[0]) / clock.record.tau0) + 1 : clock.record?.n ?? 0;

// Sliding-window deviation map for one clock (window = record / frac).
export function dynamic(clock, signal) {
  return run(request(clock, "dynamic", [
    ...rangeArgs(clock),
    "--grid", state.grid,
    "--dev", state.dyn.dev,
    "--window", Math.max(16, Math.floor(samples(clock) / state.dyn.frac)),
  ]), signal);
}

// The holdover test of one clock (training stretch = record / frac, or with
// frac 0 as many holdover times as the fit has terms: the test behind TDEV
// and HTDEV; a Kalman filter trains on a share only), with the error its
// fitted noise model predicts or the filter claims.
const holdoverRequest = (clock, frac = state.hold.frac || (state.hold.kalman ? 4 : 0)) => request(clock, "holdover", [
  ...rangeArgs(clock),
  ...(frac ? ["--train", Math.max(4, Math.floor(samples(clock) / frac)) * clock.record.tau0]
    : ["--grow", state.hold.drift ? 3 : 2]),
  ...(state.hold.drift ? ["--fit-drift"] : []),
  ...(state.hold.kalman ? ["--kalman"] : [...(state.hold.anchor ? [] : ["--no-anchor"]), "--fit"]),
  ...(state.hold.limit ? ["--limit", state.hold.limit] : []),
  "--confidence", state.hold.confidence,
]);

// Settings that decide a clock's holdover test; a change means recompute.
export const holdoverKey = (clock) => JSON.stringify(holdoverRequest(clock));

// The core's JSON as is.
export async function holdover(clock, signal) {
  const doc = await post(holdoverRequest(clock), signal);
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  return doc;
}

// The time constant the Steering page starts from: 1/128 of the record.
export const steerDefault = (clock) => Math.max(4, Math.round(samples(clock) / 128)) * clock.record.tau0;

// The steering replay of one clock: gains from a time constant or a cost, a
// steer every so many samples, the clock model fitted to the first
// record / frac alone.
const steerRequest = (clock) => {
  const o = state.steer, tau0 = clock.record.tau0;
  return request(clock, "steer", [
    ...rangeArgs(clock),
    ...(o.from === "cost" ? ["--cost", o.cost.join(":")]
      : ["--time-constant", o.T ?? steerDefault(clock), ...(o.integral ? ["--integral"] : [])]),
    ...(o.every > 1 ? ["--interval", o.every * tau0] : []),
    "--train", Math.max(16, Math.floor(samples(clock) / o.frac)) * tau0,
    ...(o.drift ? ["--fit-drift"] : []),
    "--grid", state.grid, "--series", 2000,
  ]);
};

// Settings that decide a clock's replay; a change means recompute.
export const steerKey = (clock) => JSON.stringify(steerRequest(clock));

// The core's `steer` object as is.
export async function steer(clock, signal) {
  const doc = await post(steerRequest(clock), signal);
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  return doc.steer;
}

// A timescale from several clocks; the core's JSON as is.
export async function timescale(body, signal) {
  const doc = await post(body, signal);
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  return doc;
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

// The σ(τ) that noise models predict (MDEV, ADEV and MHDEV per model), closed form.
export async function modelPreview({ n, tau0, models }, signal) {
  const args = ["model", "--synth", n, "--tau0", tau0, ...models.flatMap((m) => ["--noise-model", m.join(":")]), "--format", "json"];
  const doc = await post({ args: args.map(String), file: null }, signal);
  if (doc.schema_version !== 1) throw new Error(`the analysis core speaks schema ${doc.schema_version}; this page needs 1`);
  return doc.model;
}
