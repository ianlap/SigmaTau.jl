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
using JSON

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
               calc_ci::Bool=false,
               m_max::Union{Int, Nothing}=nothing)
    warmup()
    println("\n=== $(basename(path)) ===")
    pd, N, tau0 = load_phase_2col(path)
    ms = bench_m_values(N)
    if m_max !== nothing
        ms = ms[ms .<= m_max]
    end
    @printf("  N = %d, τ₀ = %.6g s, %d m values (%d … %d)\n",
            N, tau0, length(ms), ms[1], ms[end])
    @printf("  threads = %d, calc_ci = %s\n", Threads.nthreads(), calc_ci)
    println("  ", "-"^46)

    results = NamedTuple{(:kernel, :time, :bytes, :gctime), Tuple{Symbol, Float64, Int64, Float64}}[]
    for (sym, fn) in KERNELS
        sym in kernels || continue
        # Force a GC before timing so we don't attribute prior garbage to this kernel.
        GC.gc()
        r = @timed fn(pd, ms; calc_ci=calc_ci)
        push!(results, (kernel=sym, time=r.time, bytes=r.bytes, gctime=r.gctime))
        @printf("  %-9s  %10.3f s   alloc=%10.1f MiB   gc=%5.2fs\n",
                sym, r.time, r.bytes/2^20, r.gctime)
    end

    println("  ", "-"^46)
    @printf("  total: %.3f s\n", sum(r.time for r in results))
    return results
end

function _write_json(path::AbstractString, results, meta::Dict)
    payload = Dict(
        "meta" => meta,
        "results" => [Dict(
            "kernel" => String(r.kernel),
            "time_s" => r.time,
            "bytes" => r.bytes,
            "gctime_s" => r.gctime,
        ) for r in results],
    )
    open(path, "w") do io
        JSON.print(io, payload, 2)
    end
    println("wrote $path")
end

"Run all KERNELS on each .txt file in `synth_dir`, recording per-realization timings.
Returns a Dict suitable for JSON dump with per-kernel arrays of measurements."
function bench_synth(synth_dir::AbstractString;
                     tau0::Float64=1.0,
                     m_max::Union{Int, Nothing}=nothing,
                     kernels::Vector{Symbol}=Symbol[k for (k, _) in KERNELS],
                     calc_ci::Bool=false)
    warmup()
    files = sort([joinpath(synth_dir, f) for f in readdir(synth_dir)
                  if endswith(f, ".txt")])
    isempty(files) && error("no .txt files in $synth_dir")
    # Probe N from the first file
    pd0, N, _ = load_phase_1col(first(files); tau0=tau0)
    ms = bench_m_values(N)
    if m_max !== nothing
        ms = ms[ms .<= m_max]
    end
    println("\n=== synthetic bench ===")
    @printf("  dir = %s\n  reals = %d, N = %d, τ₀ = %.6g s\n",
            synth_dir, length(files), N, tau0)
    @printf("  %d m values (%d … %d)\n", length(ms), ms[1], ms[end])
    @printf("  threads = %d, calc_ci = %s\n", Threads.nthreads(), calc_ci)
    println("  ", "-"^46)

    # Per-kernel: vector of per-realization measurements.
    per_kernel = Dict{Symbol, Vector{NamedTuple}}()
    for (sym, _) in KERNELS
        sym in kernels || continue
        per_kernel[sym] = NamedTuple[]
    end

    for (i, path) in enumerate(files)
        pd, _, _ = load_phase_1col(path; tau0=tau0)
        for (sym, fn) in KERNELS
            sym in kernels || continue
            GC.gc()
            r = @timed fn(pd, ms; calc_ci=calc_ci)
            push!(per_kernel[sym], (
                realization = i - 1,
                time_s = r.time,
                bytes = r.bytes,
                gctime_s = r.gctime,
            ))
        end
        @printf("  realization %3d / %3d done\n", i, length(files))
    end

    println("  ", "-"^46)
    for (sym, _) in KERNELS
        sym in kernels || continue
        ts = [m.time_s for m in per_kernel[sym]]
        @printf("  %-9s  mean=%9.4f s   median=%9.4f s   std=%9.4f s\n",
                sym, sum(ts)/length(ts), sort(ts)[ceil(Int, length(ts)/2)],
                sqrt(sum((t - sum(ts)/length(ts))^2 for t in ts) / max(length(ts)-1, 1)))
    end
    return per_kernel, N, ms
end

function _write_synth_json(out::AbstractString, per_kernel, N::Int, ms::Vector{Int},
                            synth_dir::AbstractString, tau0::Float64)
    payload = Dict(
        "meta" => Dict(
            "mode" => "synth",
            "synth_dir" => synth_dir,
            "n_reals" => length(first(values(per_kernel))),
            "N" => N,
            "tau0_s" => tau0,
            "m_values" => collect(ms),
            "n_m" => length(ms),
            "m_min" => ms[1],
            "m_max" => ms[end],
            "threads" => Threads.nthreads(),
            "julia_version" => string(VERSION),
        ),
        "results" => [Dict(
            "kernel" => String(sym),
            "per_realization" => [Dict(
                "realization" => m.realization,
                "time_s" => m.time_s,
                "bytes" => m.bytes,
                "gctime_s" => m.gctime_s,
            ) for m in per_kernel[sym]],
        ) for (sym, _) in KERNELS if haskey(per_kernel, sym)],
    )
    open(out, "w") do io
        JSON.print(io, payload, 2)
    end
    println("wrote $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Two modes:
    #   <path> [out.json] [m_max]           — single-record bench (legacy)
    #   --synth <dir> [out.json] [m_max]    — synthetic statistical bench
    length(ARGS) >= 1 || error("usage:\n" *
        "  julia --project=validation/bench -t auto bench_sigmatau.jl <path> [out.json] [m_max]\n" *
        "  julia --project=validation/bench -t auto bench_sigmatau.jl --synth <dir> [out.json] [m_max]")
    if ARGS[1] == "--synth"
        length(ARGS) >= 2 || error("--synth requires a directory argument")
        synth_dir = ARGS[2]
        out_json = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, "results_sigmatau_synth.json")
        m_max_cap = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : nothing
        per_kernel, N, ms = bench_synth(synth_dir; m_max=m_max_cap)
        _write_synth_json(out_json, per_kernel, N, ms, synth_dir, 1.0)
        exit(0)
    end
    path = ARGS[1]
    out_json = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "results_sigmatau.json")
    m_max_cap = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : nothing
    pd, N, tau0 = load_phase_2col(path)
    ms = bench_m_values(N)
    if m_max_cap !== nothing
        ms = ms[ms .<= m_max_cap]
    end
    results = bench(path; m_max=m_max_cap)
    meta = Dict(
        "file" => basename(path),
        "N" => N,
        "tau0_s" => tau0,
        "n_m" => length(ms),
        "m_min" => ms[1],
        "m_max" => ms[end],
        "threads" => Threads.nthreads(),
        "julia_version" => string(VERSION),
    )
    _write_json(out_json, results, meta)
end
