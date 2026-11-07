// Estimator catalog: id, label, which axis it plots on, and a one-line note
// shown in the "+ Add estimator" picker.
export const ESTIMATORS = [
  { id: "adev",     label: "ADEV",     axis: "sigma_y", note: "Allan deviation" },
  { id: "mdev",     label: "MDEV",     axis: "sigma_y", note: "Modified Allan; separates white and flicker PM" },
  { id: "hdev",     label: "HDEV",     axis: "sigma_y", note: "Hadamard; insensitive to linear drift" },
  { id: "mhdev",    label: "MHDEV",    axis: "sigma_y", note: "Modified Hadamard" },
  { id: "pdev",     label: "PDEV",     axis: "sigma_y", note: "Parabolic deviation" },
  { id: "totdev",   label: "TOTDEV",   axis: "sigma_y", note: "Total Allan; better confidence at long τ" },
  { id: "mtotdev",  label: "MTOTDEV",  axis: "sigma_y", note: "Modified total" },
  { id: "htotdev",  label: "HTOTDEV",  axis: "sigma_y", note: "Hadamard total" },
  { id: "mhtotdev", label: "MHTOTDEV", axis: "sigma_y", note: "Modified Hadamard total" },
  { id: "theo1",    label: "Thêo1",    axis: "sigma_y", note: "Theo1; reaches τ up to 75 % of the record" },
  { id: "theoh",    label: "ThêoH",    axis: "sigma_y", note: "Hybrid of Allan and Theo1" },
  { id: "tdev",     label: "TDEV",     axis: "sigma_x", note: "Time deviation" },
  { id: "htdev",    label: "HTDEV",    axis: "sigma_x", note: "Hadamard time deviation" },
  { id: "ttotdev",  label: "TTOTDEV",  axis: "sigma_x", note: "Total time deviation" },
  { id: "mtie",     label: "MTIE",     axis: "sigma_x", note: "Maximum time interval error" },
  { id: "tierms",   label: "TIE rms",  axis: "sigma_x", note: "RMS time interval error" },
];

export const byId = Object.fromEntries(ESTIMATORS.map((e) => [e.id, e]));
export const label = (id) => byId[id]?.label ?? id.toUpperCase();

export const GRIDS = [
  ["octave", "Octave"], ["half-octave", "Half octave"], ["quarter-octave", "Quarter octave"],
  ["decade", "Decade"], ["half-decade", "Half decade"], ["all", "All τ"],
];

export const CONFIDENCES = [[0.683, "68.3 %"], [0.9, "90 %"], [0.95, "95 %"], [0.99, "99 %"]];

// SigmaTau's noise symbols (:WHPM … :RWFM) shown the usual way.
export const NOISE_NAMES = {
  WHPM: "WPM", FLPM: "FPM", WHFM: "WFM", FLFM: "FFM", RWFM: "RWFM",
};
export const noiseLabel = (s) => NOISE_NAMES[String(s)] ?? String(s);
