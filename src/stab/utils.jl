# utils.jl — Shared helpers for the stability API

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

| `kernel`                                                   | m_max          |
|------------------------------------------------------------|----------------|
| `:adev`, `:totdev`, `:pdev`                                | `(N − 1) ÷ 2`  |
| `:mdev`, `:tdev`, `:mtotdev`, `:ttotdev`, `:htotdev`       | `N ÷ 3`        |
| `:hdev`                                                    | `(N − 1) ÷ 3`  |
| `:mhdev`, `:htdev`, `:mhtotdev`                            | `N ÷ 4`        |
| `:mtie`                                                    | `N − 1`        |

Throws `ArgumentError` for unknown kernel symbols or `N` too short to
admit any `m ≥ 1`.
"""
function _default_m_values(N::Int, kernel::Symbol)
    m_max = if kernel === :adev || kernel === :totdev || kernel === :pdev
        (N - 1) ÷ 2
    elseif kernel === :mdev || kernel === :tdev || kernel === :mtotdev ||
           kernel === :ttotdev || kernel === :htotdev
        N ÷ 3
    elseif kernel === :hdev
        (N - 1) ÷ 3
    elseif kernel === :mhdev || kernel === :htdev || kernel === :mhtotdev
        N ÷ 4
    elseif kernel === :mtie
        N - 1
    else
        throw(ArgumentError("_default_m_values: unknown kernel symbol :$kernel"))
    end
    m_max < 1 && throw(ArgumentError(
        "_default_m_values: N=$N is too short to support any m for :$kernel"))
    return [2^k for k in 0:floor(Int, log2(m_max))]
end
