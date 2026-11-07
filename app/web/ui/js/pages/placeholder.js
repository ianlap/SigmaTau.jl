// Spectrum and Dynamic pages come in later milestones; the core does not
// return those yet.
const TEXT = {
  spectrum: ["Spectrum", "S<sub>y</sub>(f), S<sub>x</sub>(f) and ℒ(f), overlaid per clock."],
  dynamic: ["Dynamic deviation", "σ(t, τ) over sliding windows: one map per clock, side by side on a shared colour scale."],
};

export function placeholder(page) {
  return {
    mount(root) {
      const [title, body] = TEXT[page];
      root.innerHTML = `<div class="placeholder"><div class="card"><h2>${title}</h2><p>${body}</p><p class="small">Coming in a later step.</p></div></div>`;
    },
    refresh() {},
    plotElement() { return null; },
    hint() { return ""; },
  };
}
