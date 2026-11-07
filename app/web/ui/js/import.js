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

const isNum = (x) => x !== "" && Number.isFinite(Number(x));

export function sniff(name, text) {
  const lines = text.split(/\r?\n/);
  if (lines.length && lines[lines.length - 1] === "") lines.pop();
  const [split, delimName] = splitter(name);

  let comments = 0;
  while (comments < lines.length && lines[comments].trim().startsWith("#")) comments++;

  // Non-numeric lines after the comments are header lines the reader skips.
  let header = 0;
  while (comments + header < lines.length && header < 40) {
    const f = split(lines[comments + header]);
    if (f.length && f.every(isNum)) break;
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
    if (!bad && !f.every(isNum)) bad = j + 1;
    for (let c = 0; c < f.length; c++) (columns[c] ??= []).push(Number(f[c]));
    nRows++;
  }
  const nCols = columns.length;
  // "# tau0: 1" (SigmaTau's save) or "Tau: 1" (Stable32-style headers).
  const tauComment = lines.slice(0, start).map((c) => /\btau0?\s*[:=]\s*([0-9.eE+-]+)/i.exec(c)).find(Boolean);
  const names = nameFields.length === nCols ? nameFields : [];

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
    tau0Hint: tauComment ? Number(tauComment[1]) : null,
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

// Median step and gaps of a time column, scanning the whole file.
export function timeStats(col, unit = 1) {
  const d = [];
  for (let i = 1; i < col.length; i++) d.push(col[i] - col[i - 1]);
  if (!d.length) return null;
  const sorted = [...d].sort((a, b) => a - b);
  const step = sorted[Math.floor(sorted.length / 2)];
  let gaps = 0, missing = 0, uneven = 0;
  for (const x of d) {
    if (x > 1.5 * step) { gaps++; missing += Math.round(x / step) - 1; }
    else if (Math.abs(x - step) > 0.01 * Math.abs(step)) uneven++;
  }
  return { tau0: step * unit, gaps, missing, even: gaps === 0 && uneven === 0, span: (col[col.length - 1] - col[0]) * unit };
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
