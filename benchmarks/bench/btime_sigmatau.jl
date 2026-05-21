# benchmarks/bench/btime_sigmatau.jl
# Minimal BenchmarkTools-based bench for the three deviations most users
# reach for first: overlapping ADEV, MDEV, overlapping HDEV. Mirrors the
# Python script `btime_allantools.py` so the numbers line up directly.
#
# Synthetic randn() phase data; no file IO. Single threaded BLAS (we don't
# use BLAS in the kernels — this just keeps it from grabbing background
# cores). Total-family kernels are NOT included here; they are heavyweight
# and live in the existing `bench_sigmatau.jl` wall-clock bench.
#
# Usage (from the persistent REPL, with `benchmarks/bench/Project.toml`
# activated):
#
#     include("benchmarks/bench/btime_sigmatau.jl")
#     btime_run()                # N = 2^15
#     btime_run(N = 1 << 17)     # bigger record
#
# Or one-shot from the shell:
#
#     julia --project=benchmarks/bench benchmarks/bench/btime_sigmatau.jl
#     julia --project=benchmarks/bench benchmarks/bench/btime_sigmatau.jl 131072

using BenchmarkTools
using LinearAlgebra
using Printf
using Random
using SigmaTau

BLAS.set_num_threads(1)

const KERNELS = [
    (:adev, adev),   # overlapping ADEV  ↔ allantools oadev
    (:mdev, mdev),   # modified ADEV    ↔ allantools mdev
    (:hdev, hdev),   # overlapping HDEV  ↔ allantools ohdev
]

_m_grid(N::Integer) = [1 << k for k in 0:floor(Int, log2(N / 3))]

function btime_run(; N::Integer = 1 << 15, seed::Integer = 0)
    rng = MersenneTwister(seed)
    pd  = PhaseData(randn(rng, N), 1.0)
    ms  = _m_grid(N)

    @printf("N = %d, %d m values (%d … %d), threads = %d\n",
            N, length(ms), ms[1], ms[end], Threads.nthreads())
    println("-"^54)
    for (sym, fn) in KERNELS
        @printf("%-6s ", sym)
        @btime $fn($pd, $ms; calc_ci=false)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1 << 15
    btime_run(; N=N)
end
