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
