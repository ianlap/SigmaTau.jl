// Number and duration formatting shared by the pages.
export function sci(x, digits = 4) {
  if (x === null || x === undefined || !Number.isFinite(x)) return "–";
  if (x === 0) return "0";
  return x.toExponential(digits - 1);
}

export function num(x, digits = 4) {
  if (x === null || x === undefined || !Number.isFinite(x)) return "–";
  if (Math.abs(x) >= 1e5 || Math.abs(x) < 1e-3) return sci(x, digits);
  return Number(x.toPrecision(digits)).toString();
}

export function thousands(n) {
  return n.toLocaleString("en-US").replace(/,/g, " ");
}

export function duration(s) {
  if (!Number.isFinite(s)) return "–";
  if (s < 120) return `${num(s, 3)} s`;
  if (s < 2 * 3600) return `${num(s / 60, 3)} min`;
  if (s < 3 * 86400) return `${num(s / 3600, 3)} h`;
  return `${num(s / 86400, 3)} d`;
}

export function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
}

// A Kalman filter's clock model: flicker levels only where the fit found them.
export const modelText = (m) => [["R", m.R], ["q₁", m.q1], ["q₂", m.q2], ["h₋₁", m.flicker_fm], ["h₁", m.flicker_pm]]
  .filter(([, v], i) => i < 3 || v > 0).map(([k, v]) => `${k} = ${sci(v, 3)}`).join(", ");

// The core's refusals name command-line flags and library functions; a page says
// them in its own words. Anything else is shown as the core wrote it.
const ALGORITHM = { greenhall: "Greenhall", jones_tryon: "Jones–Tryon", stein: "Stein" };
const SAID = [
  [/^--algorithm (\w+) fits a model per clock and needs three or more clocks/, (_, a) => `${ALGORITHM[a] ?? a} fits a model to each clock and needs three or more clocks.`],
  [/^(?:timescale: )?wmax below (1\/\d+) leaves weights that cannot sum to 1/, (_, k) => `The largest weight is below ${k}, so the weights of the clocks cannot sum to 1. Raise it to ${k} or more.`],
  [/^(?:powerlaw_fit: )?(\d+) usable τ points for (\d+) components; pass fewer `alphas` or a longer record/, (_, n, k) => `The noise model has ${n} usable τ points for ${k} noise types: the record, or the stretch it is fitted to, is too short or has no readings.`],
  [/^(?:steer: )?the steady state is worked out for a steer at every sample; the gains' interval is ([\d.e+-]+) s and tau0 is ([\d.e+-]+) s/, (_, g, t) => `The steady state is worked out for a steer at every sample, and this loop steers every ${num(g / t, 3)} samples.`],
];
export function inWords(message) {
  for (const [re, say] of SAID) if (re.test(message)) return message.replace(re, say);
  return message.replace(/^(?:holdover|steer|timescale): /, "");   // the command's own prefix
}

// A UTC start time for the record check: seconds to the microsecond, the float noise of
// the core's timestamp ("…01.007999897Z") rounded away, no trailing zeros.
export function startText(iso) {
  const m = /^(\d{4}-\d\d-\d\d)T(\d\d:\d\d:\d\d)(?:\.(\d+))?(Z?)$/.exec(iso);
  if (!m) return iso.replace("T", " ");
  const us = Math.round(Number(`0.${m[3] ?? 0}`) * 1e6);
  const t = new Date(Date.parse(`${m[1]}T${m[2]}Z`) + Math.floor(us / 1e6) * 1000).toISOString().slice(0, 19).replace("T", " ");
  const frac = String(us % 1e6).padStart(6, "0").replace(/0+$/, "");
  return `${t}${frac ? `.${frac}` : ""}${m[4]}`;
}
