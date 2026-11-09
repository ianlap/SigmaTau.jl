// Estimator catalog. `family` is the short name shown in the panel rows,
// `note` the one-line description in the Add estimators dialog.
// id, label, group and axis must match `estimators()` (checked in test/stab/registry.jl).
export const GROUPS = [
  { id: "allan", title: "Allan family", sub: "σ<sub>y</sub>" },
  { id: "total", title: "Total", sub: "tighter intervals at long τ" },
  { id: "theo",  title: "Thêo", sub: "τ up to 3/4 of the record" },
  { id: "time",  title: "Time", sub: "σ<sub>x</sub> in seconds" },
];

export const ESTIMATORS = [
  { id: "adev",     label: "ADEV",     group: "allan", axis: "sigma_y", family: "Allan",                    note: "Allan, the standard" },
  { id: "mdev",     label: "MDEV",     group: "allan", axis: "sigma_y", family: "Modified Allan",           note: "separates white and flicker PM" },
  { id: "hdev",     label: "HDEV",     group: "allan", axis: "sigma_y", family: "Hadamard",                 note: "ignores linear drift" },
  { id: "mhdev",    label: "MHDEV",    group: "allan", axis: "sigma_y", family: "Modified Hadamard",        note: "modified Hadamard" },
  { id: "pdev",     label: "PDEV",     group: "allan", axis: "sigma_y", family: "Parabolic",                note: "parabolic, for averaging counters" },
  { id: "totdev",   label: "TOTDEV",   group: "total", axis: "sigma_y", family: "Total Allan",              note: "total Allan" },
  { id: "mtotdev",  label: "MTOTDEV",  group: "total", axis: "sigma_y", family: "Total modified Allan",     note: "total modified Allan" },
  { id: "htotdev",  label: "HTOTDEV",  group: "total", axis: "sigma_y", family: "Total Hadamard",           note: "total Hadamard" },
  { id: "mhtotdev", label: "MHTOTDEV", group: "total", axis: "sigma_y", family: "Total modified Hadamard",  note: "total modified Hadamard" },
  { id: "theo1",    label: "THEO1",    group: "theo",  axis: "sigma_y", family: "Thêo1",                    note: "Thêo1, bias-corrected" },
  { id: "theoh",    label: "THEOH",    group: "theo",  axis: "sigma_y", family: "ThêoH",                    note: "Thêo1 joined to ADEV" },
  { id: "tdev",     label: "TDEV",     group: "time",  axis: "sigma_x", family: "Time deviation",           note: "time deviation" },
  { id: "htdev",    label: "HTDEV",    group: "time",  axis: "sigma_x", family: "Hadamard time",            note: "Hadamard time deviation" },
  { id: "ttotdev",  label: "TTOTDEV",  group: "time",  axis: "sigma_x", family: "Total time",               note: "total time deviation" },
  { id: "mtie",     label: "MTIE",     group: "time",  axis: "sigma_x", family: "Max time interval error",  note: "max time interval error" },
  { id: "tierms",   label: "TIE rms",  group: "time",  axis: "sigma_x", family: "RMS time interval error",  note: "rms time interval error" },
];

export const byId = Object.fromEntries(ESTIMATORS.map((e) => [e.id, e]));
export const label = (id) => byId[id]?.label ?? id.toUpperCase();

export const GRIDS = [
  ["octave", "Octave"], ["half-octave", "Half octave"], ["quarter-octave", "Quarter octave"],
  ["decade", "Decade"], ["half-decade", "Half decade"], ["all", "All τ"],
];
export const gridName = (g) => (GRIDS.find(([v]) => v === g) ?? [g, g])[1];

export const CONFIDENCES = [[0.683, "68.3 %"], [0.9, "90 %"], [0.95, "95 %"], [0.99, "99 %"]];

// SigmaTau's noise symbols shown the usual way, with the power-law exponent α
// of S_y(f) ∝ f^α.
const NOISE_NAMES = { WHPM: "WPM", FLPM: "FPM", WHFM: "WFM", FLFM: "FFM", RWFM: "RWFM", FWFM: "FWFM", RRFM: "RRFM" };
export const NOISE_ALPHA = { WHPM: 2, FLPM: 1, WHFM: 0, FLFM: -1, RWFM: -2, FWFM: -3, RRFM: -4 };
export const NOISE_LONG = { WHPM: "white phase", FLPM: "flicker phase", WHFM: "white frequency", FLFM: "flicker frequency",
  RWFM: "random-walk frequency", FWFM: "flicker-walk frequency", RRFM: "random-run frequency" };
export const noiseLabel = (s) => NOISE_NAMES[String(s)] ?? String(s);
export const noiseAlpha = (s) => NOISE_ALPHA[String(s)];
const alphaText = (a) => (a > 0 ? `+${a}` : a < 0 ? `−${-a}` : "0");
export const noiseTitle = (s) => (s in NOISE_ALPHA ? `${NOISE_LONG[s]} noise, α = ${alphaText(NOISE_ALPHA[s])}` : String(s));

// Slope of log σ against log τ for a power-law noise. The modified estimators
// (MDEV, MHDEV, MTOTDEV, MHTOTDEV) separate white from flicker phase noise;
// the others give τ^−1 for both.
const MODIFIED = new Set(["mdev", "mhdev", "mtotdev", "mhtotdev"]);
export function noiseSlope(dev, noise) {
  const a = NOISE_ALPHA[String(noise)];
  if (a === undefined) return null;
  return MODIFIED.has(dev) ? -(a + 1) / 2 : -(Math.min(a, 1) + 1) / 2;
}

// Estimators the bias-correction option changes.
export const BIAS_CORRECTED = new Set(["totdev", "mtotdev", "htotdev", "mhtotdev", "ttotdev", "theo1"]);
