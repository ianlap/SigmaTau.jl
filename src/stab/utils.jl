# utils.jl — Shared helpers for the stability API

"""
    _f64(v::AbstractVector{<:Real}) → Vector{Float64}

Promote a sample vector to the `Vector{Float64}` the core kernels and
`identify_noise` require. A no-op (returns the same array, no copy) when `v`
is already a `Vector{Float64}`; otherwise allocates a converted copy. Lets
`PhaseData{Float32}` (and any other `AbstractFloat` element type) flow through
the public deviation API instead of hitting a `MethodError` on the
`Vector{Float64}`-typed kernels.
"""
_f64(v::Vector{Float64}) = v
_f64(v::AbstractVector{<:Real}) = Vector{Float64}(v)

"""
    _freq_to_phase(data::FrequencyData) → PhaseData

Convert fractional-frequency samples to phase residuals using the running
integral `x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]`. Length is preserved (N → N), matching the
convention that `adev(FrequencyData(y, τ₀), …)` agrees with
`adev(PhaseData(cumsum(y)·τ₀, τ₀), …)`.
"""
_freq_to_phase(data::FrequencyData) = PhaseData(cumsum(data.y) .* data.tau0, data.tau0)

"""
    _phase_to_freq(data::PhaseData) → FrequencyData

Convert phase residuals to fractional frequency via the canonical first
difference `y[k] = (x[k+1] − x[k]) / τ₀`. Length drops from `N` to
`N − 1`. The inverse direction (`_freq_to_phase`) recovers an
`N`-length record but loses the absolute phase offset, so the
phase ↔ frequency map is bijective only up to a constant: round-trip
`_freq_to_phase ∘ _phase_to_freq` returns `x[2:end] .- x[1]`.

Deviation estimators are shift-invariant (their second/third differences
annihilate the offset), so e.g. `adev(pd, m_values) ≈
adev(_phase_to_freq(pd), m_values)` up to the `N → N−1` length change.
"""
_phase_to_freq(data::PhaseData) = FrequencyData(diff(data.x) ./ data.tau0, data.tau0)

"""
    _default_m_values(N::Int, kernel::Symbol) → Vector{Int}

Octave-spaced default averaging-factor grid `[1, 2, 4, …, 2^k]` bounded
above by the kernel's algorithmic m-max — the largest `m` for which the
core L-check still yields at least one window. Used by the zero-arg
convenience methods on every public deviation API
(`adev(pd)`, `mdev(pd)`, …) so callers can skip the `m_values` argument.

Per-kernel m-max (derived from the `L`/`Ne` guard in each `_*_core`):

| `kernel`                  | m_max          |
|---------------------------|----------------|
| `:adev`, `:pdev`          | `(N − 2) ÷ 2`  |
| `:totdev`                 | `(N − 1) ÷ 2`  |
| `:mdev`, `:tdev`          | `(N − 1) ÷ 3`  |
| `:mtotdev`, `:ttotdev`    | `N ÷ 3`        |
| `:hdev`                   | `(N − 2) ÷ 3`  |
| `:htotdev`                | `(N − 1) ÷ 3`  |
| `:mhdev`, `:htdev`        | `(N − 1) ÷ 4`  |
| `:mhtotdev`               | `N ÷ 4`        |
| `:mtie`                   | `N − 1`        |

The ordinary and modified estimators require ≥2 analysis windows at the
largest τ (a single window is one raw difference, EDF ≈ 1), so their caps
sit one step below the total family, whose subsequence extension keeps a
single subsequence a valid long-τ estimate. (HTOTDEV's general branch runs
on `y = diff(x)` of length `N−1`, giving `(N−1) ÷ 3`; MTOTDEV runs on phase
directly, giving `N ÷ 3`.)

Throws `ArgumentError` for unknown kernel symbols or `N` too short to
admit any `m ≥ 1`.
"""
_default_m_values(N::Int, kernel::Symbol) = _grid(Octave, _kernel_m_max(N, kernel))
