// Build the 40 x 32 in, three-column poster from content.js and figures/.
// Usage: node build_poster.js   -> SigmaTau_poster_draft.pptx  (prints per-column fill)
const path = require('path');
const fs = require('fs');
const pptxgen = require('pptxgenjs');   // `npm install` in poster/, or set NODE_PATH
const C = require('./content.js');

// ---------- geometry (inches) ----------
const W = 40, H = 32, M = 1.0, MB = 0.8, GUT = 0.8, NCOL = 3;
const COLW = (W - 2 * M - (NCOL - 1) * GUT) / NCOL;   // 12.13
const PAD = 0.35;                                       // inside a panel card
const INNER = COLW - 2 * PAD;
const TITLE_H = 4.0;                                    // crimson title band
const COL_TOP = TITLE_H + 0.6;
const REF_H = 1.45;                                     // references strip
const COL_BOTTOM = H - MB - REF_H - 0.2;
const GAP = 0.45;                                       // between panels

// ---------- type scale (pt) ----------
const T = { title: 60, subtitle: 28, authors: 30, affil: 24, heading: 33, body: 24, caption: 19.5, table: 18, code: 19, refs: 15 };
const LH = 1.15;
const FONT = 'Calibri', MONO = 'Courier New';

// ---------- colors (hex, no #) ----------
const CRIMSON = '9E1B32', INK = '0B0B0B', INK2 = '52514E', RULE = 'D9D4D0', CARD = 'FFFFFF', PAGE = 'FFFFFF', TINT = 'F7F2F3', CODEBG = 'F4F4F2';

function pngSize(file) { const b = fs.readFileSync(file); return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) }; }
const plain = t => t.replace(/[\^_]\{([^}]*)\}/g, '$1');
function lines(text, widthIn, pt, em = 0.45) {
  const cpl = Math.max(6, Math.floor((widthIn * 72) / (pt * em)));
  return plain(text).split('\n').reduce((n, para) => n + Math.max(1, Math.ceil(para.length / cpl)), 0);
}
const inH = (nLines, pt) => (nLines * pt * LH) / 72;
// '_{x}' -> subscript run, '^{x}' -> superscript run
function runs(text, opts = {}) {
  return text.split(/(\^\{[^}]*\}|_\{[^}]*\})/).filter(Boolean).map(pt => {
    if (pt.startsWith('^{')) return { text: pt.slice(2, -1), options: { ...opts, superscript: true } };
    if (pt.startsWith('_{')) return { text: pt.slice(2, -1), options: { ...opts, subscript: true } };
    return { text: pt, options: { ...opts } };
  });
}
function paraRuns(text, firstOpts, lastOpts) {
  const r = runs(text); Object.assign(r[0].options, firstOpts); Object.assign(r[r.length - 1].options, lastOpts); return r;
}

const pres = new pptxgen();
pres.defineLayout({ name: 'POSTER', width: W, height: H });
pres.layout = 'POSTER';
pres.theme = { headFontFace: FONT, bodyFontFace: FONT };
pres.author = C.authors; pres.title = 'SigmaTau.jl poster draft';
const slide = pres.addSlide();
slide.background = { color: PAGE };

// ---------- title band ----------
slide.addShape(pres.ShapeType.rect, { x: 0, y: 0, w: W, h: TITLE_H, fill: { color: CRIMSON }, line: { color: CRIMSON } });
slide.addText(runs(C.title), { x: M, y: 0.35, w: W - 2 * M, h: 1.35, fontSize: T.title, bold: true, color: 'FFFFFF', fontFace: FONT, valign: 'middle', margin: 0, isTextBox: true, objectName: 'Title', fit: 'shrink' });
slide.addText(runs(C.subtitle), { x: M, y: 1.8, w: W - 2 * M, h: 1.2, fontSize: T.subtitle, color: 'FFFFFF', fontFace: FONT, valign: 'middle', margin: 0, isTextBox: true, objectName: 'Subtitle' });
slide.addText([{ text: C.authors, options: { bold: true, fontSize: T.authors } }, { text: '     ' + C.affiliation, options: { fontSize: T.affil } }],
  { x: M, y: 3.15, w: W - 2 * M, h: 0.6, color: 'FFFFFF', fontFace: FONT, valign: 'middle', margin: 0, isTextBox: true, objectName: 'Authors' });

// ---------- block renderers: each returns the height it consumed ----------
function heading(x, y, text) {
  const h = inH(lines(text, INNER, T.heading, 0.46), T.heading) + 0.05;
  slide.addText(runs(text), { x: x + PAD, y, w: INNER, h, fontSize: T.heading, bold: true, color: CRIMSON, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Heading' });
  return h + 0.18;
}
function paragraphs(x, y, paras) {
  const rr = []; let n = 0;
  paras.forEach((p, i) => { rr.push(...paraRuns(p, { paraSpaceAfter: 8 }, { breakLine: i < paras.length - 1 })); n += lines(p, INNER, T.body); });
  const h = inH(n, T.body) + (paras.length - 1) * 8 / 72 + 0.1;
  slide.addText(rr, { x: x + PAD, y, w: INNER, h, fontSize: T.body, color: INK, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Body' });
  return h + 0.12;
}
function bullets(x, y, items) {
  const rr = []; let n = 0;
  items.forEach(([lead, text], i) => {
    rr.push({ text: lead, options: { bold: true, bullet: { indent: 22 }, paraSpaceAfter: 7 } });
    rr.push(...paraRuns(text, {}, { breakLine: i < items.length - 1 }));
    n += lines(lead + text, INNER - 0.3, T.body, 0.44);
  });
  const h = inH(n, T.body) + (items.length - 1) * 7 / 72 + 0.1;
  slide.addText(rr, { x: x + PAD, y, w: INNER, h, fontSize: T.body, color: INK, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Bullets' });
  return h + 0.12;
}
function code(x, y, linesOfCode) {
  const h = linesOfCode.length * T.code * 1.3 / 72 + 0.3;
  slide.addShape(pres.ShapeType.roundRect, { x: x + PAD, y, w: INNER, h, fill: { color: CODEBG }, line: { color: RULE, width: 0.75 }, rectRadius: 0.08 });
  const rr = linesOfCode.map((l, i) => ({ text: l, options: { breakLine: i < linesOfCode.length - 1 } }));
  slide.addText(rr, { x: x + PAD + 0.2, y: y + 0.12, w: INNER - 0.4, h: h - 0.24, fontSize: T.code, fontFace: MONO, color: INK, valign: 'middle', margin: 0, isTextBox: true, objectName: 'Code' });
  return h + 0.15;
}
function figure(x, y, names, figH) {
  const n = names.length, gap = 0.3, slotW = (INNER - (n - 1) * gap) / n;
  let maxH = 0;
  names.forEach((name, i) => {
    const file = path.join(__dirname, 'figures', name + '.png');
    const { w, h } = pngSize(file);
    let fh = figH, fw = fh * w / h;
    if (fw > slotW) { fw = slotW; fh = fw * h / w; }
    slide.addImage({ path: file, x: x + PAD + i * (slotW + gap) + (slotW - fw) / 2, y, w: fw, h: fh, objectName: name });
    maxH = Math.max(maxH, fh);
  });
  return maxH + 0.12;
}
function caption(x, y, text) {
  const h = inH(lines(text, INNER, T.caption), T.caption) + 0.05;
  slide.addText(runs(text), { x: x + PAD, y, w: INNER, h, fontSize: T.caption, color: INK2, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Caption' });
  return h + 0.1;
}
function table(x, y, tb) {
  const ncol = tb.header.length;
  const base = tb.colW || (ncol === 7 ? [0.9, 2.6, 1.5, 1.6, 1.4, 1.4, 1.4] : Array(ncol).fill(1));
  const scale = INNER / base.reduce((a, b) => a + b, 0);
  const cw = base.map(c => c * scale);
  const all = [tb.header, ...tb.rows];
  const rowH = all.map(r => Math.max(T.table * 1.6 / 72, Math.max(...r.map((c, j) => lines(c, cw[j] - 0.12, T.table, 0.47))) * T.table * 1.2 / 72 + 0.12));
  const rows = [
    tb.header.map(t => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: CRIMSON }, align: 'center' } })),
    ...tb.rows.map((r, i) => r.map((t, j) => ({ text: t, options: { align: (tb.leftCols || [1]).includes(j) ? 'left' : 'center', fill: { color: i % 2 ? TINT : CARD } } }))),
  ];
  slide.addTable(rows, { x: x + PAD, y, w: INNER, colW: cw, rowH, fontSize: T.table, fontFace: FONT, color: INK, border: { type: 'solid', pt: 0.75, color: RULE }, margin: 0.05, valign: 'middle', objectName: tb.name || 'Table' });
  return rowH.reduce((a, b) => a + b, 0) + 0.15;
}

// ---------- lay the columns out ----------
const report = [];
C.columns.forEach((blocks, ci) => {
  const x = M + ci * (COLW + GUT);
  let y = COL_TOP;
  blocks.forEach((b) => {
    const y0 = y;
    const cardIdx = slide._slideObjects.length;
    slide.addShape(pres.ShapeType.rect, { x, y: y0, w: COLW, h: 1, fill: { color: CARD }, line: { color: RULE, width: 1 } });
    y += PAD;
    y += heading(x, y, b.heading);
    if (b.body) y += paragraphs(x, y, b.body);
    if (b.bullets) y += bullets(x, y, b.bullets);
    if (b.code) y += code(x, y, b.code);
    if (b.table && b.tableAfter) { if (b.tableLead) y += caption(x, y, b.tableLead) - 0.05; y += table(x, y, b.table); }
    if (b.figure) y += figure(x, y, [b.figure], b.figH);
    if (b.figures) y += figure(x, y, b.figures, b.figH);
    if (b.caption) y += caption(x, y, b.caption);
    if (b.table && !b.tableAfter) { if (b.tableLead) y += caption(x, y, b.tableLead) - 0.05; y += table(x, y, b.table); }
    y += PAD - 0.1;
    slide._slideObjects[cardIdx].options.h = y - y0;
    report.push(`col ${ci + 1}  ${b.heading.slice(0, 48).padEnd(50)} ${y0.toFixed(2)} -> ${y.toFixed(2)}`);
    y += GAP;
  });
  const end = y - GAP;
  report.push(`col ${ci + 1} ends at ${end.toFixed(2)} in  (limit ${COL_BOTTOM.toFixed(2)})  ${end > COL_BOTTOM ? 'OVERFLOW by ' + (end - COL_BOTTOM).toFixed(2) : 'ok, slack ' + (COL_BOTTOM - end).toFixed(2)}`);
});

// ---------- references strip ----------
const ry = H - MB - REF_H;
slide.addShape(pres.ShapeType.line, { x: M, y: ry - 0.15, w: W - 2 * M, h: 0, line: { color: RULE, width: 1.5 } });
const refA = C.references.slice(0, 4), refB = C.references.slice(4);
const mk = (arr, off) => arr.map((r, i) => ({ text: `${i + 1 + off}. ${r}`, options: { breakLine: i < arr.length - 1, paraSpaceAfter: 2 } }));
const refW = (W - 2 * M) * 0.36;
slide.addText([{ text: 'References', options: { bold: true, color: CRIMSON, breakLine: true } }, ...mk(refA, 0)],
  { x: M, y: ry, w: refW, h: REF_H, fontSize: T.refs, color: INK2, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'References A' });
slide.addText([{ text: ' ', options: { breakLine: true } }, ...mk(refB, 4)],
  { x: M + refW + 0.3, y: ry, w: refW, h: REF_H, fontSize: T.refs, color: INK2, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'References B' });
slide.addText(C.acknowledgment, { x: M + 2 * refW + 0.6, y: ry, w: (W - 2 * M) - 2 * refW - 0.6, h: REF_H, fontSize: T.refs, color: INK2, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Acknowledgments' });

pres.writeFile({ fileName: path.join(__dirname, 'SigmaTau_poster_draft.pptx') }).then(() => { console.log(report.join('\n')); console.log('wrote SigmaTau_poster_draft.pptx'); });
