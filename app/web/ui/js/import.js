// Client-side look at a data file for the import dialog. It follows the rules
// of SigmaTau's readers (leading # comment lines skipped, delimiter from the
// extension, else runs of whitespace) so the preview shows the columns the
// core will see; the real parse happens in the core.

const PREVIEW_ROWS = 9;

function splitter(name, lines) {
  const ext = name.toLowerCase().split(".").pop();
  let delim = ext === "csv" ? "," : ext === "tsv" ? "\t" : null;
  // A blank-separated file named .csv still reads: the extension's delimiter
  // counts only when the first data line holds it, as in the readers.
  const first = lines.find((l) => l.trim() && !l.trim().startsWith("#"));
  if (delim && first !== undefined && !first.includes(delim)) delim = null;
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
    : names.some((n) => /freq|\bffe\b|osc|fractional|ratio/i.test(n)) ? "frequency" : null;

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

// With no names to go by, readings that never exceed 1e-11 in size can only be
// fractional frequency (or a ratio): a time error that small is not a clock.
function sizeKind(columns) {
  const med = (c) => { const v = c.filter(Number.isFinite).map(Math.abs).sort((a, b) => a - b); return v.length ? v[v.length >> 1] : NaN; };
  const varies = (c) => c.length > 1 && c.some((v) => v !== c[0]);
  const readings = columns.filter((c) => varies(c) && !c.every((v, i) => i === 0 || v > c[i - 1]));
  return readings.length && readings.every((c) => med(c) < 1e-11) ? "frequency" : null;
}

// A BIPM UTC - UTC(k) file ("MJD  [UTC-UTC(NIST)]/ns  uA/ns …"): an integer MJD and the
// difference in ns per row, then uncertainty and note columns that vary by row, so the
// generic column reader cannot take it; the core reads it by its layout (bipm_utc).
const BIPM_HEAD = /^\s*MJD\s+\[UTC-UTC\((.+?)\s*\)\]\s*\/\s*ns/i;
function sniffBipm(name, lines, lab) {
  const rows = [], mjd = [], ns = [];
  for (const l of lines) {
    const f = l.trim().split(/\s+/);
    if (f.length < 2 || !/^\d+$/.test(f[0]) || !Number.isFinite(Number(f[1]))) continue;
    if (rows.length < PREVIEW_ROWS) rows.push({ ln: lines.indexOf(l) + 1, f: [f[0], f[1]] });
    mjd.push(Number(f[0])); ns.push(Number(f[1]));
  }
  const head = lines.findIndex((l) => BIPM_HEAD.test(l));
  return {
    name, lines: lines.length, delimName: "whitespace-separated", comments: 0, header: 0, namesLine: 0,
    names: ["MJD", `UTC-UTC(${lab}) (ns)`], top: [{ ln: head + 1, text: lines[head].trim() }], rows,
    nRows: mjd.length, nCols: 2, columns: [mjd, ns], bad: null, tau0Hint: null,
    format: "bipm_utc", formatLabel: `BIPM UTC − UTC(${lab})`, kindHint: "phase", stable32: false,
  };
}

// A RINEX 3 clock file (IGS / CODE clock products): "AS G08 2025 08 17 00 00 00.000000  2  bias  sigma"
// per satellite (AS) or station (AR) clock and epoch. One column per clock, on the union of the
// epochs (MJD days, NaN where a clock has no record); the core reads one clock by name (rinex_clk).
const RINEX_HEAD = /^\s*\d+\.\d+\s+C\b.*RINEX VERSION \/ TYPE/;
function sniffRinex(name, lines) {
  const body = lines.findIndex((l) => l.includes("END OF HEADER")) + 1;
  const epochs = new Map(), mjd = [], series = new Map();
  for (let j = body; j < lines.length; j++) {
    const l = lines[j], kind = l.slice(0, 3);
    if (kind !== "AS " && kind !== "AR ") continue;
    const f = l.trim().split(/\s+/);
    if (f.length < 10) continue;
    const key = f.slice(2, 8).join(" ");
    let e = epochs.get(key);
    if (e === undefined) {
      e = mjd.length; epochs.set(key, e);
      mjd.push(Date.UTC(+f[2], f[3] - 1, +f[4], +f[5], +f[6]) / 86400000 + 40587 + Number(f[7]) / 86400);
      for (const c of series.values()) c.col.push(NaN);
    }
    let c = series.get(f[1]);
    if (!c) { c = { col: new Array(mjd.length).fill(NaN), kind: kind === "AS " ? "satellite" : "station", n: 0 }; series.set(f[1], c); }
    c.col[e] = Number(f[9]); c.n++;
  }
  const names = [...series.keys()].sort();
  return {
    name, lines: lines.length, delimName: "whitespace-separated", comments: 0, header: body, namesLine: 0,
    names: ["MJD", ...names], top: [], rows: [], nRows: mjd.length, nCols: names.length + 1,
    columns: [mjd, ...names.map((n) => series.get(n).col)], bad: null, tau0Hint: null,
    format: "rinex_clk", formatLabel: "RINEX clock file", kindHint: "phase", stable32: false,
    clocks: names.map((n) => `${series.get(n).kind} · ${series.get(n).n.toLocaleString("en-US")} epochs`),
  };
}

export function sniff(name, text) {
  const lines = text.split(/\r?\n/);
  if (lines.length && lines[lines.length - 1] === "") lines.pop();
  if (lines.slice(0, 5).some((l) => RINEX_HEAD.test(l))) return sniffRinex(name, lines);
  const bipm = lines.slice(0, 5).map((l) => BIPM_HEAD.exec(l)).find(Boolean);
  if (bipm) return sniffBipm(name, lines, bipm[1]);
  const [split, delimName] = splitter(name, lines);

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
  // "# tau0: 1" (SigmaTau's save), "Tau: 1" (Stable32-style headers) or
  // "# data interval 1.0 s" (allantools' sample files).
  const tauComment = lines.slice(0, start).map((c) => /\btau0?\s*[:=]\s*([0-9.eE+-]+)/i.exec(c)
    ?? /^#\s*(?:data|sample) interval:?\s+([0-9.eE+-]+)\s*s\s*$/i.exec(c)).find(Boolean);
  // A last comment line of names ("# t  ratio  flag", as the ROCIT files have) names the columns too.
  const commentNames = comments > 0 && !header ? lines[comments - 1].replace(/^#\s*/, "").trim().split(/\s+/) : [];
  const names = timelab ? [] : nameFields.length === nCols ? nameFields
    : commentNames.length === nCols && commentNames.every((x) => !isField(x)) ? commentNames : [];
  const fmt = formatOf(name, timelab, header, rows, columns);
  const interval = timelab && head.map((l) => /"Sample Interval"\s+([0-9.eE+-]+)/.exec(l)).find(Boolean);

  return {
    name,
    lines: lines.length,
    delimName,
    comments,
    header,
    namesLine: nameFields.length === nCols && !timelab ? start : 0,       // 1-based line number of the column names
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
    kindHint: fmt?.format ? fmt.format === "table" ? tableKind(names) : fmt.kind ?? null : tableKind(names) ?? sizeKind(columns),
    stable32: /stable32/i.test(lines.slice(0, start).join("\n")) || /\.dat$/i.test(name) && header > 0,
  };
}

// The scale a column name states: "(ns)", "[ps]" for time, "ppb" for frequency.
const NAME_UNITS = [[/[(\[]\s*ns\s*[)\]]/i, 1e-9], [/[(\[]\s*(?:us|µs|μs)\s*[)\]]/i, 1e-6], [/[(\[]\s*ms\s*[)\]]/i, 1e-3],
  [/[(\[]\s*ps\s*[)\]]/i, 1e-12], [/ppb/i, 1e-9], [/ppm/i, 1e-6], [/ppt/i, 1e-12]];
const unitOf = (name) => (NAME_UNITS.find(([re]) => re.test(name)) ?? [null, 1])[1];

// Guess which columns hold time and value: a strictly increasing first column
// is time; the value is the next column that is not constant.
export function guessRoles(s) {
  if (s.clocks) return { time: 1, value: 2, timeUnit: 86400, scaling: 1 };   // MJD, then the first clock
  if (s.nCols <= 1) return { time: 0, value: 1, timeUnit: 1, scaling: unitOf(s.names[0] ?? "") };
  const inc = (c) => c.length > 1 && c.every((v, i) => i === 0 || v > c[i - 1]);
  const isMjd = (c, name) => /mjd/i.test(name) || (c[0] > 30000 && c[0] < 90000 && c[1] - c[0] < 1);
  const cand = [];
  for (let c = 0; c < s.nCols; c++) if (inc(s.columns[c])) cand.push(c + 1);
  const nm = (c) => s.names[c - 1] ?? "";
  // A clock column that only drifts upward is increasing too, so names count:
  // a column named for a phase or frequency is never time.
  const reading = (c) => /phase|freq|ratio|value|offset|clock|^[xy]$/i.test(nm(c));
  const timeLike = cand.filter((c) => !reading(c));
  // A column named like elapsed seconds wins; then the first increasing
  // column that is not MJD; then the first one, read in days when it is MJD.
  let time = timeLike.find((c) => /elapsed|sec|time|^t$|_s$/i.test(nm(c)) && !isMjd(s.columns[c - 1], nm(c)))
    ?? timeLike.find((c) => !isMjd(s.columns[c - 1], nm(c))) ?? timeLike[0] ?? 0;
  // Unnamed columns: the first one is time when it is increasing.
  if (!time && cand[0] === 1) time = 1;
  const timeUnit = time && isMjd(s.columns[time - 1], nm(time)) ? 86400 : 1;
  const named = s.names.findIndex((n, i) => i + 1 !== time && !timeLike.includes(i + 1) && /phase|freq|value|^x$|^y$/i.test(n));
  if (named >= 0) return { time, value: named + 1, timeUnit, scaling: unitOf(nm(named + 1)) };
  let value = 1;
  while (value <= s.nCols && (value === time || timeLike.includes(value))) value++;
  if (value > s.nCols) value = time === 1 ? 2 : 1;
  return { time, value, timeUnit, scaling: unitOf(nm(value)) };
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

// A typed scale: a number or a short expression with pi, + - * / and
// parentheses, such as 1/(2*pi*518.295836590863e12) (a phase in radians of a
// carrier). NaN when it is empty, malformed or not a positive number.
export function parseFactor(text) {
  const tok = text.toLowerCase().match(/\s*(\d+\.?\d*(?:e[+-]?\d+)?|\.\d+(?:e[+-]?\d+)?|pi|[-+*/()]|\S)/gy)?.map((t) => t.trim()) ?? [];
  let i = 0;
  const expr = () => {
    let v = term();
    while (tok[i] === "+" || tok[i] === "-") v = tok[i++] === "+" ? v + term() : v - term();
    return v;
  };
  const term = () => {
    let v = unary();
    while (tok[i] === "*" || tok[i] === "/") v = tok[i++] === "*" ? v * unary() : v / unary();
    return v;
  };
  const unary = () => (tok[i] === "-" ? (i++, -unary()) : tok[i] === "+" ? (i++, unary()) : atom());
  const atom = () => {
    const t = tok[i++];
    if (t === "(") { const v = expr(); return tok[i++] === ")" ? v : NaN; }
    return t === "pi" ? Math.PI : /^[\d.]/.test(t ?? "") ? Number(t) : NaN;
  };
  const v = tok.length ? expr() : NaN;
  return i === tok.length && Number.isFinite(v) && v > 0 ? v : NaN;
}

// Settings of a loaded clock hold one number for the scale. Back in the dialog
// it is a menu unit, or else a typed factor and goes in the factor field.
export function scaleFields(kind, scaling, factor = "") {
  const menu = kind === "phase" ? PHASE_UNITS : FREQ_UNITS;
  return menu.some(([f]) => f === scaling) || factor ? { scaling, factor } : { scaling: 1, factor: String(scaling) };
}

export const PHASE_UNITS = [[1, "seconds"], [1e-3, "milliseconds"], [1e-6, "microseconds"], [1e-9, "nanoseconds"], [1e-12, "picoseconds"]];
export const FREQ_UNITS = [[1, "fractional (Δf/f)"], [1e-6, "parts per million"], [1e-9, "parts per billion"], [1e-12, "parts per trillion"]];
export const TIME_UNITS = [[1, "Time (s)"], [86400, "Time (days, MJD)"], [3600, "Time (h)"], [60, "Time (min)"]];
