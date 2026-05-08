# Total-family `detrend` kwarg — design

**Date:** 2026-05-07
**Status:** Approved (brainstorm complete; ready for implementation plan)
**Track:** A4 (stability correctness)
**Branch:** `track-a4-detrend`

## Context

Cross-validation against external numerical references (Stable32, allantools)
on the four total-family kernels (`totdev`, `mtotdev`, `htotdev`, `mhtotdev`)
shows methodology-level disagreements that the existing test corpus
absorbs as wide rtol floors:

| Kernel | Stable32 rtol floor (current) | Source of disagreement |
|---|---|---|
| TOTDEV | 0.15 | Global LS detrend on top of mean-flip reflection (SigmaTau adds drift removal that SP1065 eqn 25 / Howe 1995 does not) |
| MTOTDEV | 0.05 | Bias-policy column (kernel matches) |
| HTOTDEV | 0.10 | Bias-policy column (kernel matches) |
| MHTOTDEV | n/a | No external reference exists; estimator is novel to SigmaTau |

We want to be able to reproduce each external reference's methodology
exactly, instead of tolerating a wide rtol band that hides regressions.
That means making the detrend / boundary-extension recipe a user-selectable
kwarg, with the canonical method per kernel as the default.

Pre-1.0; behavioral changes in default output are acceptable when flagged
in the changelog with migration guidance.

## Goals

1. Add a `detrend::Symbol` kwarg to all four total-family kernels (core +
   API wrappers).
2. Ship four recipes: `:howe`, `:greenhall`, `:linear`, `:legacy`.
3. Default per kernel matches the canonical methodology for that kernel.
4. Preserve bit-for-bit parity contract via opt-in `detrend=:legacy`.
5. Tighten Stable32 cross-validation TOTDEV rtol from 0.15 → 1e-4 under
   the new `:howe` default.

## API surface

Kwarg name: `detrend`. Lives on both layers:

- Core kernels: `_totdev_core`, `_mtotdev_core`, `_htotdev_core`,
  `_mhtotdev_core`.
- Public API: `totdev`, `mtotdev`, `htotdev`, `mhtotdev`.

Public API passes through to core; core does the dispatch.

### Recipe set (4 values)

| Recipe | Detrend strategy | Extension | Canonical for |
|---|---|---|---|
| `:howe` | none | endpoint mean-flip reflection: `x*_{1−j} = 2x_1 − x_{1+j}` | TOTDEV (Howe 1995, NIST SP1065 eqn 25) |
| `:greenhall` | per-window half-mean slope | per-window time-reverse: `[reverse \| original \| reverse]` | MTOTDEV / HTOTDEV (Greenhall 2003) |
| `:linear` | per-window full LS fit | per-window time-reverse | (research / current MHTOTDEV pattern) |
| `:legacy` | per-kernel current SigmaTau choice | per-kernel current SigmaTau choice | Backward-compat parity contract |

### Defaults per kernel

| Kernel | Default | Rationale |
|---|---|---|
| `totdev` | `:howe` | SP1065 eqn 25 canonical |
| `mtotdev` | `:greenhall` | Greenhall 2003 canonical |
| `htotdev` | `:greenhall` | Greenhall 2003 canonical (on frequency series) |
| `mhtotdev` | `:greenhall` | No external canonical exists; align with the Hadamard-modified family |

### Validation

Unknown recipe throws `ArgumentError("unknown detrend recipe: $recipe; valid: :howe, :greenhall, :linear, :legacy")`. Explicit, no silent fallthrough.

## Recipe semantics (math)

### `:howe` — no detrend, mean-flip reflection

Build extended series `x*` of length `3N − 4`:

```
x*_{1−j} = 2·x_1 − x_{1+j}        for j = 1..N−2     (left mean-flip)
x*_i     = x_i                     for i = 1..N       (original centered)
x*_{N+j} = 2·x_N − x_{N−j}         for j = 1..N−2     (right mean-flip)
```

Apply the kernel's standard 2nd / 3rd-difference operator on `x*`. No drift
removal anywhere. Source: Howe 1995, NIST SP1065 eqn 25; matches allantools'
default `totdev`.

### `:greenhall` — per-window half-mean slope, time-reverse reflection

For each starting index `n` and averaging factor `m`, take a `3m`-sample
window. Estimate slope from the half-mean difference:

```
half = floor(3m / 2)
slope = (mean(x[n+half : n+3m−1]) − mean(x[n : n+half−1])) / (half · τ₀)
```

Detrend the window by subtracting `slope · τ₀ · (j − 1)`. Build a
3-segment time-reverse extension of length `3·(3m)`:

```
ext = [reverse(detrended_window) | detrended_window | reverse(detrended_window)]
```

Apply the kernel operator on `ext`. Source: Greenhall 2003; matches Stable32 /
allantools for MTOT/HTOT.

For HTOTDEV the same procedure applies, but operating on the frequency
series `y = diff(x) / τ₀` rather than phase. The kernel then uses a
third-difference operator.

### `:linear` — per-window full LS detrend, time-reverse reflection

Same per-window structure as `:greenhall`, but instead of a half-mean
slope estimate, fit a full least-squares line via closed-form analytic
sums:

```
For each window x[n : n+L−1] of length L:
    y(j) = a + b·j  fit by analytic-formula LS
    detrended[j] = x[n+j−1] − (a + b·j)
```

Then the same `[reverse | original | reverse]` extension. This is the
pattern current MHTOTDEV uses internally; `:linear` extends it to the
other three kernels for research / experimental comparison.

### `:legacy` — per-kernel current behavior, asymmetric

Preserves the existing per-kernel implementations bit-for-bit:

| Kernel | Slope removal | Extension |
|---|---|---|
| TOTDEV | global LS fit over the whole vector | endpoint mean-flip reflection |
| MTOTDEV | per-window half-mean (≡ `:greenhall`) | per-window time-reverse (≡ `:greenhall`) |
| HTOTDEV | per-window half-mean on frequency (≡ `:greenhall`) | per-window time-reverse (≡ `:greenhall`) |
| MHTOTDEV | per-window full LS (≡ `:linear`) | per-window time-reverse (≡ `:linear`) |

So `:legacy` reduces to existing recipes for MTOT/HTOT/MHTOT. Only TOTDEV's
`:legacy` is genuinely distinct (it adds the global pre-detrend on top of
`:howe`).

### Default-change impact

| Kernel | Old output (current) | New default output | Change |
|---|---|---|---|
| TOTDEV | global LS + mean-flip reflect | reflect-only (Howe) | **changes** |
| MTOTDEV | half-mean + time-reverse | half-mean + time-reverse | identical |
| HTOTDEV | half-mean + time-reverse (freq) | half-mean + time-reverse (freq) | identical |
| MHTOTDEV | full LS + time-reverse | half-mean + time-reverse | **changes** |

## EDF and bias implications

The detrend recipe affects estimator-to-true-variance ratio (bias) and
degrees of freedom (EDF). Calibration status under the new defaults:

- **`bias_correction(:totvar, ...)`** uses Howe 1995 / SP1065 formula
  `B = 1 − a·(τ/T)`. Calibrated for `:howe` methodology. **Switching
  default to `:howe` fixes a long-standing miscalibration** where the
  current `:legacy` global LS detrend was already removing drift before
  the bias factor was applied.
- **`bias_correction(:mtot, ...)`** SP1065 table calibrated for
  `:greenhall`. Already correct under the new default.
- **`bias_correction(:htot, ...)`** SP1065 table calibrated for
  `:greenhall`. Already correct under the new default.
- **`bias_correction(:mhtot, ...)`** identity (no canonical model exists,
  per PR `8947445`). New default `:greenhall` doesn't change that.
- **`_coeff_mhtot` (EDF table)** empirical SP1065 fit, originally tuned
  for `:linear` MHTOTDEV. Under new `:greenhall` default the table
  becomes approximate. Caveat in docstring; full re-fit deferred to
  future Monte Carlo work.

For non-default recipes (e.g. `mtotdev(...; detrend=:linear)`) bias and
EDF were calibrated for the canonical recipe of that kernel. Docstring
caveat: "Bias and EDF are calibrated for the canonical recipe of this
kernel; non-canonical recipes are user-beware."

## Testing strategy

### Updated tests

| Test | Update | Rtol |
|---|---|---|
| `legacy_kernels.jl` parity (52 assertions) | Add `detrend=:legacy` kwarg to every `_*_core` call. | 1e-12 (unchanged) |
| Stable32 TOTDEV rows | Tighten under new `:howe` default. | 0.15 → 1e-4 |
| Stable32 MTOTDEV / HTOTDEV rows | No change. Existing 0.05 / 0.10 floors are bias-policy-driven, not kernel-driven; unchanged under this PR. | 0.05 / 0.10 (unchanged) |
| allantools cross-val TOTDEV | **Defer to follow-up PR after `track-b1-allantools` merges to main.** Same tightening as Stable32. | 0.15 → 1e-4 (in follow-up) |

### New tests (added in this PR)

1. **Cross-recipe equivalence** — `:legacy` MTOT ≡ `:greenhall` MTOT and
   `:legacy` HTOT ≡ `:greenhall` HTOT at rtol=1e-12 on the
   `legacy_kernels.jl` fixture. ~16 assertions (8 m-values × 2 kernels).
   Locks the alias claim; catches silent drift.
2. **Per-recipe smoke** — for each of `(:howe, :greenhall, :linear, :legacy)`
   × each of `(:totdev, :mtotdev, :htotdev, :mhtotdev)`, run on
   `_gen_powerlaw_phase` for one mid-spectrum noise type (α=0, WHFM).
   Assert per (recipe, kernel): all output values finite, all positive
   (no boundary-failure zeros), and within an order of magnitude of the
   `:legacy` reference output across the m-grid (catches gross
   implementation errors without overconstraining the recipe).
   ~16 assertions plus magnitude bounds.
3. **TOTDEV `:howe` ↔ Stable32 tight match** — explicit assertion that
   `:howe` TOTDEV matches Stable32's TOTDEV at rtol=1e-4 (replaces the
   old rtol=0.15 floor). Confirms the `:howe` recipe matches SP1065 eqn 25.
4. **MHTOTDEV `:greenhall` finite-output smoke** — across all 5
   power-law noise types (`_gen_powerlaw_phase` α ∈ {2, 1, 0, -1, -2}),
   assert finite, reasonable values. New default needs at least basic
   coverage.

## Migration / breaking changes

Pre-1.0; defaults change for TOTDEV and MHTOTDEV. Documented in CHANGELOG
under `[Unreleased] → Changed` with migration guidance:

> TOTDEV and MHTOTDEV default outputs change vs the previous release.
> TOTDEV now defaults to the SP1065 eqn 25 canonical methodology
> (`detrend=:howe`); previous behavior is available via
> `totdev(...; detrend=:legacy)`. MHTOTDEV now defaults to the
> Greenhall-style per-window half-mean detrend (`detrend=:greenhall`);
> previous behavior is `mhtotdev(...; detrend=:legacy)`. MTOTDEV and
> HTOTDEV default outputs are unchanged.

## Out-of-scope / future work (track in TODO.md)

- **MHTOTDEV bias / EDF Monte Carlo.** Synthesize known-noise via
  `_gen_powerlaw_phase` for each α; compute MHTOT and MHDEV; the ratio
  yields bias factor B(α). Fit `_coeff_mhtot` empirically per detrend
  recipe (`:greenhall`, `:linear`). Track per-recipe to be honest about
  EDF being recipe-specific.
- **allantools cross-val tightening for TOTDEV.** After
  `track-b1-allantools` merges to main, update
  `lib/SigmaTauStability/test/allantools_cross_validation.jl` to drop the
  TOTDEV rtol from 0.15 → 1e-4 under new `:howe` default.
- **Bias / EDF coverage for non-canonical recipe combinations.** When a
  user calls e.g. `mtotdev(...; detrend=:linear)` we currently fall back
  to the canonical-recipe SP1065 tables, with a docstring caveat. A
  separate Monte Carlo campaign could fit per-(kernel, recipe) bias and
  EDF tables.

## Implementation organization

Per-recipe internal functions (rather than one mega-function with
branches per kernel). Each recipe is a tight 30–50 lines of clear math.
A single top-level dispatcher per kernel selects which to call. Aliases
(`:legacy` MTOT, etc.) just delegate to the canonical implementation:

```
function _totdev_core(x, m_values, tau0; detrend=:howe)
    detrend === :howe     && return _totdev_howe(x, m_values, tau0)
    detrend === :greenhall && return _totdev_greenhall(x, m_values, tau0)
    detrend === :linear   && return _totdev_linear(x, m_values, tau0)
    detrend === :legacy   && return _totdev_legacy(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end
```

Legacy dispatchers for MTOT/HTOT/MHTOT delegate to the canonical recipe
implementation.

## References

- Howe, D. A. (1995). *The Total Deviation Approach to Long-Term
  Characterization of Frequency Stability.* IEEE Trans. UFFC 42(2).
  Defines TOTVAR and the mean-flip reflection extension.
- Greenhall, C. A., & Riley, W. J. (2003). *Uncertainty of Stability
  Variances Based on Finite Differences.* PTTI proceedings. Defines
  per-window half-mean detrending and EDF approximations.
- NIST Special Publication 1065 (Riley & Howe). *Handbook of Frequency
  Stability Analysis.* Authoritative reference for the deviation family
  and bias / EDF tables.
- IEEE Std 1139-2022. *Standard Definitions of Physical Quantities for
  Fundamental Frequency and Time Metrology — Random Instabilities.*
- allantools (Anders Wallin), `github.com/aewallin/allantools`. Reference
  Python implementation; `totdev` follows SP1065 eqn 25 verbatim.
