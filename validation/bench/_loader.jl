# validation/bench/_loader.jl
# Tiny phase-record loader for the 2-column .txt files under
# `reference/clock_data/`. Format: each row is `<MJD-time> <phase-value>`,
# whitespace-separated. We infer τ₀ from the median Δt of the time column
# (in seconds) and return a `PhaseData`.

using SigmaTau: PhaseData
using DelimitedFiles: readdlm
using Statistics: median

function load_phase_2col(path::AbstractString)
    arr = readdlm(path, Float64)  # Matrix{Float64} of shape (N, 2)
    size(arr, 2) == 2 || error("$(path): expected 2 columns, got $(size(arr, 2))")
    N = size(arr, 1)
    # MJD time (days) → step in seconds.
    tau0 = median(diff(@view arr[:, 1])) * 86400.0
    return PhaseData(arr[:, 2], tau0), N, tau0
end

"Per-octave m grid covering the record (caps at N/3 for total-family sanity)."
function bench_m_values(N::Integer)
    cap = floor(Int, log2(N / 3))
    return [1 << k for k in 0:cap]
end

"Single-column phase loader (used by the synthetic bench)."
function load_phase_1col(path::AbstractString; tau0::Float64=1.0)
    arr = readdlm(path, Float64)
    arr = vec(arr)
    return PhaseData(arr, tau0), length(arr), tau0
end
