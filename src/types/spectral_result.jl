"""
    SpectralResult

Unified return type for every spectral-density estimate (`Sy`, `Sx`, `L`).

Mirrors [`StabilityResult`](@ref): a flat, non-parametric record carrying the
one-sided frequency grid, the estimated spectrum, and the Welch parameters that
produced it, so a saved or plotted result is self-describing.

$(TYPEDFIELDS)
"""
struct SpectralResult
    "Which estimate produced this result: `:Sy`, `:Sx`, or `:L`."
    spectral_type::Symbol
    "One-sided frequency grid `f` in Hz."
    freq::Vector{Float64}
    "Spectral density at each `f`; interpret per `units`."
    psd::Vector{Float64}
    "Units of `psd`: `:per_Hz` (S_y, 1/Hz), `:s2_per_Hz` (S_x, s²/Hz), or `:dBc_per_Hz` (ℒ(f))."
    units::Symbol
    "Welch segment length (samples)."
    nperseg::Int
    "Welch segment overlap (samples)."
    noverlap::Int
    "Window applied to each segment: `:hann`, `:hamming`, or `:rectangular`."
    window::Symbol
end

function Base.show(io::IO, r::SpectralResult)
    n = length(r.freq)
    rng = n == 0 ? "" : ", f∈[$(r.freq[1]), $(r.freq[end])] Hz"
    print(io, "SpectralResult(:", r.spectral_type, ", ", n, " bins", rng,
              ", ", r.units, ", nperseg=", r.nperseg, ")")
end
