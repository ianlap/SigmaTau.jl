#!/usr/bin/env julia
# Time SigmaTau.mtotdev and SigmaTau.htotdev on the shared datasets.
# JULIA_NUM_THREADS controls whether this is the single-thread headline
# number or the multi-thread bonus.

using SigmaTau
using DelimitedFiles
using Printf

const TAU0 = 1.0
const NS = (1000, 4000, 16000)

octave_ms(N; frac=1/3) = begin
    out = Int[]
    m = 1
    while m <= floor(Int, N * frac)
        push!(out, m)
        m *= 2
    end
    out
end

function time_call(f; repeats=3)
    best = Inf
    result = nothing
    for _ in 1:repeats
        t0 = time_ns()
        result = f()
        dt = (time_ns() - t0) / 1e9
        if dt < best
            best = dt
        end
    end
    return best, result
end

function main()
    here = @__DIR__
    data_dir = joinpath(here, "data")
    results_dir = joinpath(here, "results")
    isdir(results_dir) || mkpath(results_dir)
    nthreads = Threads.nthreads()
    suffix = nthreads == 1 ? "single" : "threaded$(nthreads)"
    csv_path = joinpath(results_dir, "julia_results_$(suffix).csv")
    devs_path = joinpath(results_dir, "julia_devs_$(suffix).txt")

    # warmup on a tiny dataset to charge JIT off the benchmark window
    println("# JULIA_NUM_THREADS=$(nthreads)")
    println("# warming up SigmaTau kernels...")
    warm_x = cumsum(randn(64))
    let pd = PhaseData(warm_x, TAU0)
        mtotdev(pd, [1, 2]; calc_ci=false, correct_bias=false)
        htotdev(pd, [1, 2]; calc_ci=false, correct_bias=false)
    end
    println("# warmup done")

    rows = String["N,kernel,impl,seconds"]
    devs_lines = String[]
    for N in NS
        path = joinpath(data_dir, "phase_N$(N).txt")
        x = vec(readdlm(path, Float64))
        m_values = octave_ms(N)
        taus = Float64.(m_values) .* TAU0
        pd = PhaseData(x, TAU0)
        println("\n## N=$N  m_values=$m_values")

        for (kname, fn) in (("mtotdev", mtotdev), ("htotdev", htotdev))
            label = "SigmaTau.$kname ($suffix)"
            print("  timing $label ... ")
            t, res = time_call(
                () -> fn(pd, m_values; calc_ci=false, correct_bias=false);
                repeats=3,
            )
            @printf("%.4fs\n", t)
            push!(rows, "$N,$kname,sigmatau_$(suffix),$(t)")
            push!(devs_lines, "$kname.N$N=" * join(res.dev, ","))
        end
    end

    open(csv_path, "w") do io
        for r in rows; println(io, r); end
    end
    open(devs_path, "w") do io
        for l in devs_lines; println(io, l); end
    end
    println("\nwrote $csv_path")
    println("wrote $devs_path")
end

main()
