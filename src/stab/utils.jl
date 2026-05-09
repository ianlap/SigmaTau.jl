# utils.jl — Shared helpers for the stability API

"""
    _freq_to_phase(data::FrequencyData) → PhaseData

Convert fractional-frequency samples to phase residuals using the running
integral `x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]`. Length is preserved (N → N), matching the
convention that `adev(FrequencyData(y, τ₀), …)` agrees with
`adev(PhaseData(cumsum(y)·τ₀, τ₀), …)`.
"""
_freq_to_phase(data::FrequencyData) = PhaseData(cumsum(data.y) .* data.tau0, data.tau0)
