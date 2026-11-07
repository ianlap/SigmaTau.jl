// Spectrum and Dynamic pages are next; the core does not return those yet.
const TEXT = {
  spectrum: ["Spectrum", "S<sub>y</sub>(f), S<sub>x</sub>(f) and ℒ(f) per clock. Needs a spectrum command in the core, which the FFTW trim fix unblocks."],
  dynamic: ["Dynamic deviation", "σ(t, τ) over sliding windows, as a heat map per clock. Needs a dynamic command in the core."],
};

export function mountPlaceholder(root, page) {
  const [title, body] = TEXT[page];
  root.innerHTML = `<div class="placeholder"><div class="panel empty-card"><h2>${title}</h2><p class="muted">${body}</p><p class="muted small">Coming in a later step.</p></div></div>`;
  return { refresh() {} };
}
