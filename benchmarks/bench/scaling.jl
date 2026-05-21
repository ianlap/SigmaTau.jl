# benchmarks/bench/scaling.jl
#
# Empirical scaling fits per kernel: T(N) ≈ a · N^b. Useful for sizing real
# runs ("how long will totdev take on my 10 M-sample record?") without
# actually launching them.
#
# Each kernel is benched on synthetic randn input across a range of N using
# the FULL public API (calc_ci=true, correct_bias=true defaults). The
# m-grid is the kernel's natural octave default (`_default_m_values`), so
# both the per-m work AND the number of m's grow with N — matching the
# wall-time a user actually pays for a typical call.
#
# Fit is plain log-log linear regression: log T = log a + b · log N. The
# `b` exponent reads cleanly:
#     b ≈ 1.0   → O(N) or O(N log N) — adev/mdev/totdev family
#     b ≈ 2.0   → O(N²)              — mtotdev / htotdev / mhtotdev (per-window)
#
# Emits `scaling_sigmatau.json` for combined cross-library rendering against
# allantools (see `scaling_allantools.py` + `render_scaling.py`).
#
# Usage from the REPL (with `Pkg.activate("benchmarks/bench")` first):
#
#     include("benchmarks/bench/scaling.jl")
#     fit_scaling()                       # benches + prints table + writes JSON
#     predict(:totdev,   10_000_000)      # → predicted seconds
#     predict(:mtotdev, 1_000_000)        # → predicted seconds
#
# Or one-shot from the shell:
#
#     julia --project=benchmarks/bench -t auto benchmarks/bench/scaling.jl

using BenchmarkTools
using JSON
using LinearAlgebra
using Printf
using Random
using SigmaTau

BLAS.set_num_threads(1)

# Per-kernel callable. We measure the public API with its defaults: CI on,
# bias correction on. That includes `identify_noise + calculate_edf +
# confidence_intervals` overhead — the full wall-time a user pays for a
# typical call. (allantools doesn't auto-compute CI/bias; we accept that
# asymmetry in the cross-library comparison and note it in the renderer.)
const KERNELS = [
    (:adev,     adev),
    (:mdev,     mdev),
    (:tdev,     tdev),
    (:hdev,     hdev),
    (:mhdev,    mhdev),
    (:htdev,    htdev),
    (:totdev,   totdev),
    (:mtotdev,  mtotdev),
    (:htotdev,  htotdev),
    (:mhtotdev, mhtotdev),
]

# Single N grid for the display table. The O(N²) total-family is capped
# at `SLOW_MAX` because going further takes many seconds per call.
const SLOW_KERNELS = Set([:mtotdev, :htotdev, :mhtotdev])
const NS = [1 << k for k in 10:2:18]   # 1 024, 4 096, 16 k, 65 k, 262 k
const SLOW_MAX = 65_536                # cap for SLOW_KERNELS

# Module-level cache so `predict()` works without passing fits around.
const _FITS = Ref{Dict{Symbol, NamedTuple}}(Dict{Symbol, NamedTuple}())

_humanN(N::Integer) = N >= 1024 ? "$(N >>> 10)k" : "$N"

function _fit_loglog(Ns::Vector{Int}, ts::Vector{Float64})
    lx = log.(Float64.(Ns))
    ly = log.(ts)
    n  = length(Ns)
    mx = sum(lx) / n
    my = sum(ly) / n
    Sxx = sum((lx .- mx).^2)
    Sxy = sum((lx .- mx) .* (ly .- my))
    b = Sxy / Sxx
    log_a = my - b * mx
    a = exp(log_a)
    pred = log_a .+ b .* lx
    ss_res = sum((ly .- pred).^2)
    ss_tot = sum((ly .- my).^2)
    r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 1.0
    return a, b, r2
end

"""
    fit_scaling(; seconds=1.0, predict_at=1_000_000, json_out=nothing)

Benchmark every kernel across the N grids, fit `T ≈ a · N^b` in log-log
space, print a table, cache the fits for `predict()`, and (if `json_out`
is set) write a machine-readable record for cross-library rendering.

`seconds` is the per-measurement `@belapsed` budget. Total runtime is
roughly `seconds · sum(N grid length over kernels)` plus a fixed warmup.
"""
function fit_scaling(; seconds::Float64=1.0,
                       predict_at::Int=1_000_000,
                       json_out::Union{Nothing,AbstractString}=
                           joinpath(@__DIR__, "scaling_sigmatau.json"))
    @printf("Julia threads = %d. m_values = octave default per kernel.\n",
            Threads.nthreads())
    @printf("Full API defaults (calc_ci=true, correct_bias=true).\n")
    @printf("Per-measurement budget = %.1fs.\n\n", seconds)

    # `@belapsed` parses `seconds=…` at macro-expansion time, so we set the
    # global default once instead of interpolating per call.
    prev_seconds = BenchmarkTools.DEFAULT_PARAMETERS.seconds
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = seconds

    # Header
    @printf("%-9s ", "kernel")
    for N in NS
        @printf("%10s ", "N=" * _humanN(N))
    end
    @printf("│ %-22s %5s │ %s\n",
            "T ≈ a · N^b", "R²", "predict T(N=$(_humanN(predict_at)))")
    sep_len = 11 + 11 * length(NS) + 35 + 26
    println("─"^sep_len)

    fits = Dict{Symbol, NamedTuple}()
    try
        for (sym, fn) in KERNELS
            measured_Ns = Int[]
            ts = Float64[]
            @printf("%-9s ", sym)
            for N in NS
                if sym in SLOW_KERNELS && N > SLOW_MAX
                    @printf("%10s ", "—")
                    continue
                end
                x  = randn(MersenneTwister(0), N)
                pd = PhaseData(x, 1.0)
                t  = @belapsed $fn($pd)
                push!(ts, t)
                push!(measured_Ns, N)
                @printf("%9.2es ", t)
            end
            a, b, r2 = _fit_loglog(measured_Ns, ts)
            T_at = a * Float64(predict_at)^b
            @printf("│ %.3e · N^%-5.3f  %5.3f │ %9.3f s\n", a, b, r2, T_at)
            fits[sym] = (a=a, b=b, r2=r2, Ns=measured_Ns, times=ts)
        end
    finally
        BenchmarkTools.DEFAULT_PARAMETERS.seconds = prev_seconds
    end

    _FITS[] = fits

    if json_out !== nothing
        payload = Dict(
            "meta" => Dict(
                "library" => "sigmatau",
                "julia_version" => string(VERSION),
                "threads" => Threads.nthreads(),
                "kwargs" => "calc_ci=true, correct_bias=true (defaults)",
                "predict_at" => predict_at,
                "seconds_budget" => seconds,
            ),
            "fits" => [Dict(
                "kernel" => String(sym),
                "a" => fits[sym].a,
                "b" => fits[sym].b,
                "r2" => fits[sym].r2,
                "Ns" => fits[sym].Ns,
                "times" => fits[sym].times,
            ) for (sym, _) in KERNELS if haskey(fits, sym)],
        )
        open(json_out, "w") do io
            JSON.print(io, payload, 2)
        end
        println("\nwrote $json_out")
    end

    return fits
end

"""
    predict(kernel::Symbol, N::Integer) → Float64

Predicted runtime in seconds for `kernel` on a length-`N` record, using
fits cached by the most recent `fit_scaling()` call. Extrapolation
outside the fitted N range is the user's responsibility — exponents are
robust, prefactors drift if asymptotic constants haven't settled.
"""
function predict(kernel::Symbol, N::Integer)
    isempty(_FITS[]) && error("Call fit_scaling() first to populate fits.")
    haskey(_FITS[], kernel) || error("No fit for :$kernel. Known: $(sort(collect(keys(_FITS[]))))")
    f = _FITS[][kernel]
    return f.a * Float64(N)^f.b
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_scaling()
end
