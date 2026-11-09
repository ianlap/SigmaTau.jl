// Help page: how to use the explorer, what each estimator and noise type
// means, and the keyboard shortcuts. Static text.
import { GROUPS, ESTIMATORS } from "../estimators.js";

const SECTIONS = [
  ["start", "Getting started"],
  ["pages", "The pages"],
  ["clocks", "Several clocks and a reference"],
  ["estimators", "Estimators"],
  ["noise", "Noise types"],
  ["intervals", "Confidence intervals"],
  ["drift", "Drift removal"],
  ["spectrum", "Spectra"],
  ["dynamic", "Dynamic deviation"],
  ["files", "Files"],
  ["keys", "Keyboard shortcuts"],
];

const NOISE = [
  ["WPM", "White phase modulation", "α = 2", "μ = −2", "ADEV ∝ τ⁻¹; MDEV separates it from flicker PM (MDEV ∝ τ⁻³ᐟ²)."],
  ["FPM", "Flicker phase modulation", "α = 1", "μ = −2", "ADEV ∝ τ⁻¹ like white PM; MDEV ∝ τ⁻¹."],
  ["WFM", "White frequency modulation", "α = 0", "μ = −1", "ADEV ∝ τ⁻¹ᐟ². The usual floor of a passive atomic standard."],
  ["FFM", "Flicker frequency modulation", "α = −1", "μ = 0", "ADEV flat: the flicker floor."],
  ["RWFM", "Random-walk frequency modulation", "α = −2", "μ = 1", "ADEV ∝ τ¹ᐟ². Long-term wander, often from the environment."],
];

function estimatorTable() {
  return GROUPS.map((g) => `
    <h3>${g.title} <span class="muted">· ${g.sub}</span></h3>
    <dl class="gloss">${ESTIMATORS.filter((e) => e.group === g.id).map((e) =>
      `<dt>${e.label}</dt><dd>${e.family}: ${e.note}.</dd>`).join("")}</dl>`).join("");
}

export function mount(root) {
  root.innerHTML = `
    <nav class="help-nav card" aria-label="Help sections">${SECTIONS.map(([id, t]) => `<a href="#help" data-goto="${id}">${t}</a>`).join("")}</nav>
    <article class="help card">
      <h1>Using the SigmaTau explorer</h1>

      <section id="help-start"><h2>Getting started</h2>
        <ol>
          <li>Open a phase or frequency record on <b>Data</b>, or drop a file anywhere in the window. The import dialog asks which column is time and which columns are clocks.</li>
          <li>Check the record: gaps, points past 5σ and drift are listed next to the phase plot.</li>
          <li>Pick estimators in the panel (<b>+ Add estimator</b>) and press <b>Run analysis</b> (F5).</li>
          <li>Read the result on <b>Stability</b>, <b>Spectrum</b> and <b>Dynamic</b>. Export the plot or every number from the <b>Export</b> menu.</li>
        </ol>
        <p>Nothing recomputes behind your back on Stability and Spectrum: after you change the record, the estimators or the options, the Run button and the status bar turn amber until you run again.</p>
      </section>

      <section id="help-pages"><h2>The pages</h2>
        <dl class="gloss">
          <dt>Data</dt><dd>Load a file or synthesize a test record; the record check (samples, span, gaps, outliers, drift); drift removal, gap filling and outlier removal; the phase plot.</dd>
          <dt>Stability</dt><dd>σ(τ) on log–log axes for σ<sub>y</sub> (frequency) or σ<sub>x</sub> (time) estimators, with confidence intervals as bars or a band (<b>CI band</b>) and the identified noise types along the top (<b>Noise ID</b>). One clock shows the per-τ table; several clocks show "Compare at τ" with the most stable clock in bold. Click the plot to move τ.</dd>
          <dt>Spectrum</dt><dd>S<sub>y</sub>(f) and either ℒ(f) at a carrier you choose or S<sub>x</sub>(f), on one frequency axis. Settings here apply at once.</dd>
          <dt>Dynamic</dt><dd>σ(t, τ) in sliding windows, to see whether the stability changes over the record. Recomputes when you change the estimator or window.</dd>
        </dl>
      </section>

      <section id="help-clocks"><h2>Several clocks and a reference</h2>
        <p>Add clocks with <b>+ Add clock</b>, or mark several columns as Clock in the import dialog. With more than one clock, colour means clock (eight colours, then they repeat) and line style means estimator. The eye button hides a clock from the plots.</p>
        <p>Clocks read from one file usually share the instrument's reference ("Measured against" in the import dialog). <b>Reference</b> in the Clocks list re-references them to one of the clocks, x<sub>i</sub> − x<sub>R</sub>; the reference itself leaves the plots. Clocks from other files have no shared reference or time stamps and stay as measured.</p>
      </section>

      <section id="help-estimators"><h2>Estimators</h2>
        <p>Every estimator works on phase x(t) (frequency records are integrated first) at averaging times τ = m·τ₀ from the τ grid in Options. σ<sub>y</sub> estimators are dimensionless; σ<sub>x</sub> estimators are in seconds.</p>
        ${estimatorTable()}
        <p class="muted">Thêo1 and ThêoH reach τ up to three quarters of the record but take seconds on long records. The total estimators and Thêo1 are bias-corrected when "Bias-correct totals and Thêo1" is on.</p>
      </section>

      <section id="help-noise"><h2>Noise types</h2>
        <p>The noise column of the per-τ table is SigmaTau's identification at each τ (lag-1 autocorrelation where there are enough averaged samples, the B1 and R(n) ratios at long τ). It sets the equivalent degrees of freedom and so the confidence interval. S<sub>y</sub>(f) ∝ f<sup>α</sup>; σ<sub>y</sub>²(τ) ∝ τ<sup>μ</sup>.</p>
        <p><b>Noise ID</b> on Stability shows the same identification as a strip above the plot, one run per noise type (one strip per clock when several are shown). It comes from the first σ<sub>y</sub> estimator in the panel. With one clock, dashed guides give the slope each noise type implies for that estimator: τ<sup>−1</sup> for white and flicker PM on ADEV (τ<sup>−3/2</sup> and τ<sup>−1</sup> on MDEV), τ<sup>−1/2</sup> for white FM, flat for flicker FM and τ<sup>+1/2</sup> for random-walk FM. A curve that bends away from its guide is a sign of a mixed or misidentified region. Noise identification runs with the confidence intervals, so both need them on.</p>
        <table class="help-tbl"><thead><tr><th>Label</th><th>Name</th><th>S<sub>y</sub></th><th>σ²(τ)</th><th>Signature</th></tr></thead>
          <tbody>${NOISE.map((n) => `<tr><td><span class="badge">${n[0]}</span></td><td>${n[1]}</td><td>${n[2]}</td><td>${n[3]}</td><td>${n[4]}</td></tr>`).join("")}</tbody></table>
      </section>

      <section id="help-intervals"><h2>Confidence intervals</h2>
        <p>Intervals use the χ² distribution with the equivalent degrees of freedom (EDF) for the identified noise type and the number of samples at each τ, at the confidence level in Options (68.3 % is ±1σ). EDF falls as τ grows, so the last points carry wide intervals. n<sub>eff</sub> and the interval bounds are in the row tooltips and in every export.</p>
      </section>

      <section id="help-drift"><h2>Drift removal</h2>
        <p><b>Frequency drift</b> on the Data plot fits a model to the fractional frequency by least squares and reports the drift rate at the end of the record per day, with its standard error: <b>Linear</b> a + b·t, <b>Quadratic</b> a + b·t + c·t², <b>Log</b> a·ln(1 + b·t) + c (the MIL-PRF-55310 aging model) and <b>Diffusion</b> a + b·√(t + c). <b>Residuals</b> plots what the model leaves. The standard error assumes white residuals, so it is optimistic when flicker or random-walk FM dominates. Hover the rate for the fitted coefficients.</p>
        <p><b>Drift removal</b> in the record check fits the phase by least squares, over the samples that are not gaps, and removes the fit before every analysis. <b>Linear</b> removes a frequency offset (a straight line in phase); <b>Quadratic</b> also removes a linear frequency drift (a parabola in phase), which SP1065 asks for before the total estimators. The check lists what was removed, and the phase plot shows the residual. HDEV and its relatives ignore a linear drift anyway; ADEV, MDEV and the totals rise as τ<sup>+1</sup> when one is left in.</p>
      </section>

      <section id="help-spectrum"><h2>Spectra</h2>
        <p>S<sub>y</sub>(f) (1/Hz) and S<sub>x</sub>(f) (s²/Hz) are one-sided Welch estimates with a Hann window and 50 % overlap; the segment is the largest power of two up to a quarter of the record. ℒ(f) = ½·(2π·f<sub>carrier</sub>)²·S<sub>x</sub>(f), in dBc/Hz. "Log-binned" averages ten bins per decade; "Welch, every bin" shows the raw estimate. Slope guides name stretches whose slope is close to a whole power law. A record with gaps must be filled first.</p>
      </section>

      <section id="help-dynamic"><h2>Dynamic deviation</h2>
        <p>DADEV and DHDEV compute ADEV or HDEV inside a window that slides across the record by half its length. Rows are window centre times, columns are τ, and darker means less stable. A stationary clock gives the same colours down every column; a drift, step or disturbed stretch changes them. DHDEV ignores linear drift inside each window. With several clocks, each gets its own map on one shared colour scale, and <b>Compare at τ</b> below them plots every clock's deviation over time at one τ (pick it there or click a map column), with each clock's whole-record value dotted.</p>
      </section>

      <section id="help-files"><h2>Files</h2>
        <p>Text files with one record per line: .dat, .txt, .csv (commas), .tsv (tabs), anything else split on spaces. Lines starting with # are comments wherever they are; Stable32 headers ending in "# Header End" are skipped and their Tau: line gives τ₀. Without a time column you must give τ₀. NaN marks a missing sample. Instrument files are recognised by their layout: TimeLab .tim logs (τ₀ from the header), Stable32 .phd and .frd, and counter readings such as the Keysight 53230A writes (readings in Hz become fractional frequency; give the gate time as τ₀). Settings are remembered for files with the same extension and header; Import settings… on the Data page reopens the dialog.</p>
      </section>

      <section id="help-keys"><h2>Keyboard shortcuts</h2>
        <ul class="keys">
          <li><kbd>F5</kbd> or <kbd>Ctrl</kbd>+<kbd>Enter</kbd> run the analysis</li>
          <li><kbd>1</kbd> Data · <kbd>2</kbd> Stability · <kbd>3</kbd> Spectrum · <kbd>4</kbd> Dynamic · <kbd>?</kbd> Help</li>
          <li>On a plot: hover for values, drag to zoom, double-click to reset</li>
          <li><kbd>Esc</kbd> closes a dialog</li>
        </ul>
      </section>
    </article>`;
  root.addEventListener("click", (e) => {
    const id = e.target.closest("[data-goto]")?.dataset.goto;
    if (!id) return;
    e.preventDefault();
    root.querySelector(`#help-${id}`)?.scrollIntoView({ behavior: "smooth", block: "start" });
  });
}

export function refresh() {}
export function plotElement() { return null; }
export function hint() { return "Methods, glossary and shortcuts"; }
