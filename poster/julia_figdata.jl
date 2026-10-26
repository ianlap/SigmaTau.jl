# Produce poster figure data with the real package (run from the repo root):
#   julia --project=. poster/julia_figdata.jl
# Writes poster/data/suite_fixture.csv  (every deviation on the validation fixture, with noise type, EDF, CI)
#        poster/data/dynamic_map.csv    (dadev map of a synthesized record with a mid-record noise-level change)
#        poster/data/streaming_check.csv (streamed vs batch ADEV/MDEV/HDEV/MHDEV at several sample counts)
using SigmaTau, Random, DelimitedFiles, Printf
root = normpath(joinpath(@__DIR__, ".."))
out  = joinpath(@__DIR__, "data"); mkpath(out)

# 1) suite on the shipped fixture -------------------------------------------
p = read_phase(joinpath(root, "test/fixtures/validation/stable32gen.DAT"); time_col=0, value_col=1, tau0=1.0, header=10)
devs = (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev, :totdev, :mtotdev, :htotdev, :mhtotdev, :theoh)
s = stability(p; devs=devs)
open(joinpath(out, "suite_fixture.csv"), "w") do io
    println(io, "dev,tau,sigma,lo,hi,noise_type,edf,neff")
    for r in s.results, i in eachindex(r.tau)
        isnan(r.dev[i]) && continue
        lo = isempty(r.ci) ? NaN : r.ci[i].lo; hi = isempty(r.ci) ? NaN : r.ci[i].hi
        nt = isempty(r.noise_type) ? "" : String(r.noise_type[i]); edf = isempty(r.edf) ? NaN : r.edf[i]
        @printf(io, "%s,%.17g,%.17g,%.17g,%.17g,%s,%.17g,%d\n", r.deviation_type, r.tau[i], r.dev[i], lo, hi, nt, edf, r.neff[i])
    end
end
println("suite: ", join(string.(d.deviation_type) for d in s.results), ", "))

# 2) dynamic map: white FM at 1e-12 for 4096 s, then 3e-12 plus flicker FM for 4096 s, tau0 = 1 s --
rng = Xoshiro(20261026)
a = noise_gen(PhaseData, 4096, 1.0; sigma1=Dict(0 => 1e-12), rng=rng)
b = noise_gen(PhaseData, 4096, 1.0; sigma1=Dict(0 => 3e-12, -1 => 1e-12), rng=rng)
xa = a.x; xb = b.x .- b.x[1] .+ xa[end]
pd = PhaseData(vcat(xa, xb[2:end]), 1.0)
m = [1, 2, 4, 8, 16, 32, 64, 128]
d = dadev(pd, m; window=1024, step=128)
open(joinpath(out, "dynamic_map.csv"), "w") do io
    println(io, "t," * join(("tau=" * string(t) for t in d.tau), ","))
    for (k, t) in enumerate(d.t)
        println(io, string(t) * "," * join((@sprintf("%.10g", v) for v in d.dev[k, :]), ","))
    end
end
println("dynamic map: ", size(d.dev), " windows x taus, window=", d.window)

# 3) streaming vs batch exactness ---------------------------------------------
q = noise_gen(PhaseData, 6000, 1.0; sigma1=Dict(2 => 1e-12, -2 => 2e-14), rng=Xoshiro(7))
mv = [1, 3, 10, 30, 100]
open(joinpath(out, "streaming_check.csv"), "w") do io
    println(io, "dev,n,m,streamed,batch,relerr")
    for dv in (:adev, :mdev, :hdev, :mhdev)
        acc = StreamingStability(dv, 1.0, mv)
        checkpoints = Set([500, 1000, 2000, 4000, 6000])
        for (k, x) in enumerate(q.x)
            push!(acc, x)
            if k in checkpoints
                st = snapshot(acc)
                bt = getfield(SigmaTau, dv)(PhaseData(q.x[1:k], 1.0), mv; ci=false)
                for j in eachindex(mv)
                    rel = abs(st.dev[j] - bt.dev[j]) / bt.dev[j]
                    @printf(io, "%s,%d,%d,%.17g,%.17g,%.3g\n", dv, k, mv[j], st.dev[j], bt.dev[j], rel)
                end
            end
        end
    end
end
println("streaming check written")
