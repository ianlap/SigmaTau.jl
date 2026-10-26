// Build the 40 x 32 in, three-column poster draft from content.js and figures/.
// Usage: node build_poster.js   -> SigmaTau_poster_draft.pptx
const path = require('path');
const pptxgen = require('/tmp/claude-0/-home-user-SigmaTau-jl/50547bef-133e-52b8-914a-595eb70e9a4a/scratchpad/node/node_modules/pptxgenjs');
const fs = require('fs');
const C = require('./content.js');

// ---------- geometry (inches) ----------
const W = 40, H = 32, M = 1.0, GUT = 0.8, NCOL = 3;
const COLW = (W - 2 * M - (NCOL - 1) * GUT) / NCOL;   // 12.13
const PAD = 0.35;                                       // inside a panel
const INNER = COLW - 2 * PAD;
const TITLE_H = 4.0;                                    // crimson title band height
const COL_TOP = TITLE_H + 0.6;
const REF_H = 1.75;                                      // bottom strip
const MB = 0.8;                                         // bottom margin
const COL_BOTTOM = H - MB - REF_H - 0.2;
const GAP = 0.45;                                       // between panels

// ---------- type scale (pt) ----------
const T = { title: 60, subtitle: 30, authors: 30, affil: 24, heading: 33, body: 24, caption: 19.5, table: 18, refs: 15 };
const LH = 1.15;                                        // line-height factor
const FONT = 'Calibri';

// ---------- colors (hex, no #) ----------
const CRIMSON = '9E1B32', INK = '0B0B0B', INK2 = '52514E', RULE = 'D9D4D0', CARD = 'FFFFFF', PAGE = 'FFFFFF', TINT = 'F7F2F3';

// PNG dimensions from the IHDR chunk (no dependency)
function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}
// rough text metrics: Calibri averages ~0.47 em per character in running text
function lines(text, widthIn, pt, em = 0.45) {
  const cpl = Math.max(10, Math.floor((widthIn * 72) / (pt * em)));
  return text.split('\n').reduce((n, para) => n + Math.max(1, Math.ceil(para.length / cpl)), 0);
}
const inH = (nLines, pt) => (nLines * pt * LH) / 72;

const pres = new pptxgen();
pres.defineLayout({ name: 'POSTER', width: W, height: H });
pres.layout = 'POSTER';
pres.theme = { headFontFace: FONT, bodyFontFace: FONT };
pres.author = C.authors; pres.title = 'SigmaTau.jl poster draft';
const slide = pres.addSlide();
slide.background = { color: PAGE };

// ---------- title band ----------
slide.addShape(pres.ShapeType.rect, { x: 0, y: 0, w: W, h: TITLE_H, fill: { color: CRIMSON }, line: { color: CRIMSON } });
slide.addText(C.title, { x: M, y: 0.35, w: W - 2 * M, h: 2.1, fontSize: T.title, bold: true, color: 'FFFFFF', fontFace: FONT,
  valign: 'middle', margin: 0, isTextBox: true, objectName: 'Title', fit: 'shrink' });
slide.addText(C.subtitle, { x: M, y: 2.5, w: W - 2 * M, h: 0.55, fontSize: T.subtitle, color: 'FFFFFF', fontFace: FONT, italic: true,
  valign: 'middle', margin: 0, isTextBox: true, objectName: 'Subtitle' });
slide.addText([{ text: C.authors, options: { bold: true, fontSize: T.authors } }, { text: '     ' + C.affiliation, options: { fontSize: T.affil } }],
  { x: M, y: 3.15, w: W - 2 * M, h: 0.6, color: 'FFFFFF', fontFace: FONT, valign: 'middle', margin: 0, isTextBox: true, objectName: 'Authors' });

// ---------- panels ----------
const report = [];
function heading(x, y, text) {
  const n = lines(text, INNER, T.heading, 0.46);
  const h = inH(n, T.heading) + 0.05;
  slide.addText(text, { x: x + PAD, y, w: INNER, h, fontSize: T.heading, bold: true, color: CRIMSON, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Heading' });
  return h + 0.18;
}
function paragraphs(x, y, paras) {
  const runs = [];
  let n = 0;
  paras.forEach((p, i) => { runs.push({ text: p, options: { breakLine: i < paras.length - 1, paraSpaceAfter: 8 } }); n += lines(p, INNER, T.body); });
  const h = inH(n, T.body) + (paras.length - 1) * 8 / 72 + 0.1;
  slide.addText(runs, { x: x + PAD, y, w: INNER, h, fontSize: T.body, color: INK, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Body' });
  return h + 0.12;
}
function bullets(x, y, items) {
  const runs = []; let n = 0;
  items.forEach(([lead, text], i) => {
    runs.push({ text: lead, options: { bold: true, bullet: { indent: 22 }, paraSpaceAfter: 7 } });
    runs.push({ text, options: { breakLine: i < items.length - 1 } });
    n += lines(lead + text, INNER - 0.3, T.body, 0.44);
  });
  const h = inH(n, T.body) + (items.length - 1) * 7 / 72 + 0.1;
  slide.addText(runs, { x: x + PAD, y, w: INNER, h, fontSize: T.body, color: INK, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Bullets' });
  return h + 0.12;
}
function figure(x, y, names, figH) {
  const n = names.length, gap = 0.3;
  const slotW = (INNER - (n - 1) * gap) / n;
  let maxH = 0;
  names.forEach((name, i) => {
    const file = path.join(__dirname, 'figures', name + '.png');
    const { w, h } = pngSize(file);
    let fh = figH, fw = fh * w / h;
    if (fw > slotW) { fw = slotW; fh = fw * h / w; }
    const fx = x + PAD + i * (slotW + gap) + (slotW - fw) / 2;
    slide.addImage({ path: file, x: fx, y, w: fw, h: fh, objectName: name });
    maxH = Math.max(maxH, fh);
  });
  return maxH + 0.12;
}
function caption(x, y, text) {
  const n = lines(text, INNER, T.caption);
  const h = inH(n, T.caption) + 0.05;
  slide.addText(text, { x: x + PAD, y, w: INNER, h, fontSize: T.caption, color: INK2, fontFace: FONT, valign: 'top', margin: 0, isTextBox: true, objectName: 'Caption' });
  return h + 0.1;
}
function table(x, y, block) {
  const colW = [0.9, 2.6, 1.5, 1.6, 1.4, 1.4, 1.4];
  const scale = INNER / colW.reduce((a, b) => a + b, 0);
  const cw = colW.map(c => c * scale);
  const hdr = block.header.map(t => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: CRIMSON }, align: 'center' } }));
  const rows = [hdr, ...block.rows.map((r, i) => r.map((t, j) => ({ text: t, options: { align: j === 1 ? 'left' : 'center', fill: { color: i % 2 ? TINT : CARD } } })))];
  const rowH = (T.table * 1.65) / 72;   // LibreOffice/PowerPoint pad rows; size the card for the rendered height
  slide.addTable(rows, { x: x + PAD, y, w: INNER, colW: cw, rowH, fontSize: T.table, fontFace: FONT, color: INK, border: { type: 'solid', pt: 0.75, color: RULE }, margin: 0.05, valign: 'middle', objectName: 'Coefficient table' });
  return rowH * rows.length + 0.12;
}

C.columns.forEach((blocks, ci) => {
  const x = M + ci * (COLW + GUT);
  let y = COL_TOP;
  blocks.forEach((b) => {
    const y0 = y;
    // card first (drawn under the content); height patched after measuring via a second pass
    const cardIdx = slide._slideObjects.length;
    slide.addShape(pres.ShapeType.rect, { x, y: y0, w: COLW, h: 1, fill: { color: CARD }, line: { color: RULE, width: 1 } });
    y += PAD;
    y += heading(x, y, b.heading);
    if (b.body) y += paragraphs(x, y, b.body);
    if (b.bullets) y += bullets(x, y, b.bullets);
    if (b.figure) y += figure(x, y, [b.figure], b.figH);
    if (b.figures) y += figure(x, y, b.figures, b.figH);
    if (b.caption) y += caption(x, y, b.caption);
    if (b.table) { if (b.tableLead) y += caption(x, y, b.tableLead) - 0.05; y += table(x, y, b.table); }
    y += PAD - 0.1;
    slide._slideObjects[cardIdx].options.h = y - y0;
    report.push(`col ${ci + 1}  ${b.heading.slice(0, 48).padEnd(50)} ${y0.toFixed(2)} -> ${y.toFixed(2)}`);
    y += GAP;
  });
  report.push(`col ${ci + 1} ends at ${(y - GAP).toFixed(2)} in  (limit ${COL_BOTTOM.toFixed(2)})  ${y - GAP > COL_BOTTOM ? 'OVERFLOW by ' + (y - GAP - COL_BOTTOM).toFixed(2) : 'ok, slack ' + (COL_BOTTOM - y + GAP).toFixed(2)}`);
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

pres.writeFile({ fileName: path.join(__dirname, 'SigmaTau_poster_draft.pptx') }).then(() => {
  console.log(report.join('\n'));
  console.log('wrote SigmaTau_poster_draft.pptx');
});
