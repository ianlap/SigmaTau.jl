# validation/bench/bench_sigmatau.jl
# Wall-clock benchmark for every SigmaTau public deviation on a real
# clock recording. The total-family kernels are multithreaded across the
# m-loop (see `Threads.@threads :dynamic` in core/total.jl); start Julia
# with `--threads auto` (or `-t N`) to see the parallel speedup. BLAS is
# pinned to a single thread because we don't use linalg in the kernels —
# this just keeps BLAS from grabbing cores in the background.
#
# Strategy: warm-start each kernel on a tiny (2048-sample) PhaseData
# first so the JIT/precompile cost is paid out of band. Then time the
# real run with @elapsed. Reported in seconds.
#
# Run from the persistent REPL (assuming the validation/plots project
# is active — it has SigmaTauStability dev'd):
#
#     bench("reference/clock_data/6krb25apr.txt")
#     bench("reference/clock_data/6kocxounsteered.txt")
#     bench("reference/clock_data/6k27febunsteered.txt")
#
# Each call reuses the warm-started kernels, so back-to-back invocations
# don't re-pay the JIT cost.
#
# Optional: `bench(path; kernels=[:adev, :mtotdev], calc_ci=false)` to
# run a subset, e.g. when you only want MTOTDEV on the long file.

using LinearAlgebra
using SigmaTau
using Printf

include(joinpath(@__DIR__, "_loader.jl"))

BLAS.set_num_threads(1)  # we don't use BLAS in the kernels; keep it out of the way.

if Threads.nthreads() == 1
    @warn "Julia started with -t 1 (single-threaded). Total-family kernels are " *
          "Threads.@threads-parallel; restart with `julia --project=validation/plots -t auto` " *
          "to see the multithreaded speedup."
else
    @info "Julia threads: $(Threads.nthreads())"
end

const KERNELS = [
    (:adev,     adev),
    (:mdev,     mdev),
    (:tdev,     tdev),
    (:hdev,     hdev),
    (:mhdev,    mhdev),
    (:ldev,     ldev),
    (:totdev,   totdev),
    (:mtotdev,  mtotdev),
    (:htotdev,  htotdev),
    (:mhtotdev, mhtotdev),
]

const _WARMED = Ref(false)

"Run every kernel once on a 2048-sample array so the JIT pays its cost
out of band. Idempotent — repeated calls are no-ops."
function warmup()
    _WARMED[] && return nothing
    println("Warming up (2048-sample dummy run; one-time JIT cost) …")
    pd = PhaseData(randn(2048), 1.0)
    ms_warm = [1, 2, 4, 8]
    for (sym, fn) in KERNELS
        t = @elapsed fn(pd, ms_warm; calc_ci=false)
        @printf("  warm %-9s %.3fs\n", sym, t)
    end
    _WARMED[] = true
    return nothing
end

function bench(path::AbstractString;
               kernels::Vector{Symbol}=Symbol[k for (k, _) in KERNELS],
               calc_ci::Bool=false)
    warmup()
    println("\n=== $(basename(path)) ===")
    pd, N, tau0 = load_phase_2col(path)
    ms = bench_m_values(N)
    @printf("  N = %d, τ₀ = %.6g s, %d m values (1 … %d)\n",
            N, tau0, length(ms), ms[end])
    @printf("  threads = %d, calc_ci = %s\n", Threads.nthreads(), calc_ci)
    println("  ", "-"^46)

    results = Tuple{Symbol, Float64}[]
    for (sym, fn) in KERNELS
        sym in kernels || continue
        t = @elapsed fn(pd, ms; calc_ci=calc_ci)
        push!(results, (sym, t))
        @printf("  %-9s  %10.3f s\n", sym, t)
    end

    println("  ", "-"^46)
    @printf("  total: %.3f s\n", sum(t for (_, t) in results))
    return results
end
