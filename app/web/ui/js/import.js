// Client-side look at a data file for the import dialog. It follows the rules
// of SigmaTau's readers (leading # comment lines skipped, delimiter from the
// extension, else runs of whitespace) so the preview shows the columns the
// core will see; the real parse happens in the core.

const PREVIEW_ROWS = 9;

function splitter(name) {
  const ext = name.toLowerCase().split(".").pop();
  const delim = ext === "csv" ? "," : ext === "tsv" ? "\t" : null;
  const fn = (line) => {
    const t = line.trim();
    if (!t) return [];
    return delim ? t.split(delim).map((s) => s.trim().replace(/^"(.*)"$/, "$1")) : t.split(/\s+/);
  };
  return [fn, delim === "," ? "comma-separated" : delim === "\t" ? "tab-separated" : "whitespace-separated"];
}

const TIMELAB_FIELD = /^(STR|DBL|S32|BLN|S64|U32) 0x[0-9A-Fa-f]+ "/;
// Keysight / Agilent SCPI real format, as counters log readings.
const SCPI_REAL = /^[+-]\d\.\d+E[+-]\d+$/;

// The instrument layout of a file, or null for plain columns.
function formatOf(name, timelab, header, rows, columns) {
  // Named columns with ISO 8601 timestamps (GPSDO logs) need the core's
  // table reader; the plain column reader takes numbers only.
  if (header && rows.some(({ f }) => f.some((x) => ISO_TIME.test(x)))) return { format: "table", label: "Named-column table", kind: null };
  if (timelab) return { format: "timelab", label: "TimeLab file", kind: "phase" };
  const ext = name.toLowerCase().split(".").pop();
  if (ext === "phd") return { format: "stable32", label: "Stable32 phase file", kind: "phase" };
  if (ext === "frd") return { format: "stable32", label: "Stable32 frequency file", kind: "frequency" };
  if (header || !rows.length || columns.length > 2 || !rows.every(({ f }) => SCPI_REAL.test(f[f.length - 1]))) return null;
  const v = columns[columns.length - 1].filter((x) => Number.isFinite(x) && x !== 9.91e37).map(Math.abs).sort((a, b) => a - b);
  const freq = v.length && v[v.length >> 1] >= 1;
  return freq ? { format: "counter", label: "Counter readings in Hz", kind: "frequency" }
              : { format: "counter", label: "Counter time-interval readings", kind: "phase" };
}

// The record type a table's column names suggest, as the core's table reader
// picks its default column: phase names first, then frequency.
const tableKind = (names) =>
  names.some((n) => /phase|time ?err|\btie\b|pps|time ?offset|time ?diff/i.test(n)) ? "phase"
    : names.some((n) => /freq|\bffe\b|osc|fractional/i.test(n)) ? "frequency" : null;

// NaN marks a missing sample; Julia's reader parses it as a Float64.
const isNum = (x) => x !== "" && (Number.isFinite(Number(x)) || /^[+-]?nan$/i.test(x));

// ISO 8601 timestamps, as SigmaTau's _iso_seconds reads them.
const ISO_TIME = /^(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2}(?:\.\d*)?)(Z|([+-])(\d{2}):?(\d{2}))?$/;
const isField = (x) => isNum(x) || ISO_TIME.test(x);
// A field as a number: UTC seconds for a timestamp.
function fieldValue(x) {
  const m = ISO_TIME.exec(x);
  if (!m) return Number(x);
  const off = m[8] ? (m[8] === "+" ? 1 : -1) * (Number(m[9]) * 3600 + Number(m[10]) * 60) : 0;
  return Date.UTC(+m[1], m[2] - 1, +m[3], +m[4], +m[5]) / 1000 + Number(m[6]) - off;
}

export function sniff(name, text) {
  const lines = text.split(/\r?\n/);
  if (lines.length && lines[lines.length - 1] === "") lines.pop();
  const [split, delimName] = splitter(name);

  let comments = 0;
  while (comments < lines.length && lines[comments].trim().startsWith("#")) comments++;

  // Instrument formats the core reads by layout (see SigmaTau.detect_format):
  // a TimeLab log's samples start after its "TIC n" line.
  const head = lines.slice(0, 400);
  const tic = head.findIndex((l) => /^TIC \d+/.test(l));
  const timelab = tic >= 0 && head.some((l) => TIMELAB_FIELD.test(l));

  // Non-numeric lines after the comments are header lines the reader skips.
  let header = timelab ? tic + 1 - comments : 0;
  while (!timelab && comments + header < lines.length && header < 40) {
    const f = split(lines[comments + header]);
    if (f.length && f.every(isField)) break;
    header++;
  }
  const start = comments + header;
  const nameFields = header > 0 ? split(lines[start - 1]) : [];

  const columns = [];          // full numeric columns, for gap scan and preview
  const rows = [];             // first rows as text, for the table
  let nRows = 0;
  let bad = null;
  for (let j = start; j < lines.length; j++) {
    const f = split(lines[j]);
    if (!f.length) continue;
    if (rows.length < PREVIEW_ROWS) rows.push({ ln: j + 1, f });
    if (!bad && !f.every(isField)) bad = j + 1;
    for (let c = 0; c < f.length; c++) (columns[c] ??= []).push(fieldValue(f[c]));
    nRows++;
  }
  const nCols = columns.length;
  // "# tau0: 1" (SigmaTau's save) or "Tau: 1" (Stable32-style headers).
  const tauComment = lines.slice(0, start).map((c) => /\btau0?\s*[:=]\s*([0-9.eE+-]+)/i.exec(c)).find(Boolean);
  const names = nameFields.length === nCols && !timelab ? nameFields : [];
  const fmt = formatOf(name, timelab, header, rows, columns);
  const interval = timelab && head.map((l) => /"Sample Interval"\s+([0-9.eE+-]+)/.exec(l)).find(Boolean);

  return {
    name,
    lines: lines.length,
    delimName,
    comments,
    header,
    namesLine: names.length ? start : 0,       // 1-based line number of the column names
    names,
    top: lines.slice(0, start).map((t, i) => ({ ln: i + 1, text: t })).slice(-4),
    rows,
    nRows,
    nCols,
    columns,
    bad,
    tau0Hint: interval ? Number(interval[1]) : tauComment ? Number(tauComment[1]) : null,
    format: fmt?.format ?? null,      // passed to the core as --input-format
    formatLabel: fmt?.label ?? null,
    kindHint: fmt?.format === "table" ? tableKind(names) : fmt?.kind ?? null,
    stable32: /stable32/i.test(lines.slice(0, start).join("\n")) || /\.dat$/i.test(name) && header > 0,
  };
}

// Guess which columns hold time and value: a strictly increasing first column
// is time; the value is the next column that is not constant.
export function guessRoles(s) {
  if (s.nCols <= 1) return { time: 0, value: 1, timeUnit: 1 };
  const inc = (c) => c.length > 1 && c.every((v, i) => i === 0 || v > c[i - 1]);
  const isMjd = (c, name) => /mjd/i.test(name) || (c[0] > 30000 && c[0] < 90000 && c[1] - c[0] < 1);
  const cand = [];
  for (let c = 0; c < s.nCols; c++) if (inc(s.columns[c])) cand.push(c + 1);
  const nm = (c) => s.names[c - 1] ?? "";
  // A column named like elapsed seconds wins; then any increasing column
  // that is not MJD; then MJD, read in days.
  let time = cand.find((c) => /elapsed|sec|time|^t$|_s$/i.test(nm(c)) && !isMjd(s.columns[c - 1], nm(c)))
    ?? cand.find((c) => !isMjd(s.columns[c - 1], nm(c))) ?? cand[0] ?? 0;
  const timeUnit = time && isMjd(s.columns[time - 1], nm(time)) ? 86400 : 1;
  const named = s.names.findIndex((n, i) => i + 1 !== time && !cand.includes(i + 1) && /phase|freq|value|^x$|^y$/i.test(n));
  if (named >= 0) return { time, value: named + 1, timeUnit };
  let value = 1;
  while (value <= s.nCols && (value === time || cand.includes(value))) value++;
  if (value > s.nCols) value = time === 1 ? 2 : 1;
  return { time, value, timeUnit };
}

// Grid step and gaps of a time column, scanning the whole file: the same
// step as the reader's _time_grid (the mean of the regular steps, rounded to
// the shortest value that keeps every timestamp within a quarter step of its
// slot), so coarse timestamps such as MJD to six decimals give the true τ0.
export function timeStats(col, unit = 1) {
  const t = unit === 1 ? col : col.map((x) => x * unit);
  const d = [];
  for (let i = 1; i < t.length; i++) d.push(t[i] - t[i - 1]);
  if (!d.length) return null;
  const step = gridStep(t, d);
  let gaps = 0, missing = 0, uneven = 0;
  for (const x of d) {
    if (x > 1.5 * step) { gaps++; missing += Math.round(x / step) - 1; }
    else if (Math.abs(x - step) > 0.01 * Math.abs(step)) uneven++;
  }
  return { tau0: step, gaps, missing, even: gaps === 0 && uneven === 0, span: t[t.length - 1] - t[0] };
}

function gridStep(t, d) {
  let lo = Infinity;
  for (const x of d) if (x < lo) lo = x;
  let sum = 0, n = 0;
  for (const x of d) if (x < 1.5 * lo) { sum += x; n++; }
  const est = sum / n;
  if (!(lo > 0)) return est;
  const slots = (step) => {
    const s = [0];
    for (const x of d) {
      const k = x / step;
      if (Math.abs(k - Math.round(k)) > 0.25) return null;
      s.push(s[s.length - 1] + Math.round(k));
    }
    return s;
  };
  for (let digits = 1; digits <= 6; digits++) {
    const step = Number(est.toPrecision(digits));
    if (Math.abs(step - est) >= 0.01 * est) continue;
    const s = slots(step);
    if (s && s.every((k, i) => Math.abs(t[i] - t[0] - k * step) <= 0.25 * step)) return step;
  }
  const s = slots(est);
  return s && s[s.length - 1] > 0 ? (t[t.length - 1] - t[0]) / s[s.length - 1] : est;
}

// Evenly thinned copy of a column, for the small preview plot.
export function thin(col, n = 400) {
  if (col.length <= n) return col.map((v, i) => [i, v]);
  const step = col.length / n;
  const out = [];
  for (let i = 0; i < n; i++) { const k = Math.floor(i * step); out.push([k, col[k]]); }
  return out;
}

export const PHASE_UNITS = [[1, "seconds"], [1e-3, "milliseconds"], [1e-6, "microseconds"], [1e-9, "nanoseconds"], [1e-12, "picoseconds"]];
export const FREQ_UNITS = [[1, "fractional (Δf/f)"], [1e-6, "parts per million"], [1e-9, "parts per billion"], [1e-12, "parts per trillion"]];
export const TIME_UNITS = [[1, "Time (s)"], [86400, "Time (days, MJD)"], [3600, "Time (h)"], [60, "Time (min)"]];
