// Client-side preview of a data file for the import dialog. It mirrors the
// rules of SigmaTau's readers (leading # comments skipped, delimiter from the
// extension or whitespace) closely enough to show columns; the real parse
// happens in the core with the columns the user picks.

export function sniff(name, text) {
  const lines = text.split(/\r?\n/);
  const comments = [];
  let i = 0;
  while (i < lines.length && lines[i].trim().startsWith("#")) comments.push(lines[i++]);
  const ext = name.toLowerCase().split(".").pop();
  const delim = ext === "csv" ? "," : ext === "tsv" ? "\t" : null;
  const split = (line) => {
    const t = line.trim();
    if (!t) return [];
    return delim ? t.split(delim).map((s) => s.trim().replace(/^"(.*)"$/, "$1")) : t.split(/\s+/);
  };
  const body = lines.slice(i);
  // Non-numeric rows at the top are header lines the reader must skip.
  let header = 0;
  while (header < body.length && header < 20) {
    const f = split(body[header]);
    if (f.length && f.every((x) => Number.isFinite(Number(x)))) break;
    header++;
  }
  const headerNames = header > 0 ? split(body[header - 1]) : [];
  const rows = [];
  let nRows = 0;
  for (let j = header; j < body.length; j++) {
    const f = split(body[j]);
    if (!f.length) continue;
    nRows++;
    if (rows.length < 12) rows.push(f);
  }
  const nCols = Math.max(0, ...rows.map((r) => r.length));
  const tauFromComment = comments.map((c) => /tau0\s*:\s*([0-9.eE+-]+)/.exec(c)).find(Boolean);
  return {
    comments: comments.length,
    header,
    headerNames,
    rows,
    nRows,
    nCols,
    tau0Hint: tauFromComment ? Number(tauFromComment[1]) : null,
  };
}

// Median spacing of the first rows of a time column, as a τ0 hint.
export function tau0FromColumn(rows, col) {
  const t = rows.map((r) => Number(r[col - 1])).filter(Number.isFinite);
  const d = [];
  for (let i = 1; i < t.length; i++) d.push(t[i] - t[i - 1]);
  if (!d.length) return null;
  d.sort((a, b) => a - b);
  return d[Math.floor(d.length / 2)];
}
