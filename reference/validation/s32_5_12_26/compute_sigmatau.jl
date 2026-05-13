# compute_sigmatau.jl — Run SigmaTau's deviations on the s32_5_12 phase record
# and dump (deviation, AF, tau, dev, ci_lower, ci_upper, alpha) to CSV.
#
# Run from the repo root:
#     julia --project=. reference/validation/s32_5_12_26/compute_sigmatau.jl
#
# AF ladder matches Stable32's run in all_deviations.csv.
# Bias-correction settings (after the 2026-05-12 σ/√B fix):
#   - totdev:  correct_bias=true  (Stable32 applies √B unbias for TOTVAR)
#   - mtotdev: correct_bias=false (Stable32 reports raw biased MTOT)
#   - htotdev: correct_bias=true  (Stable32 empirically matches √B unbias)
# ttotdev is derived from mtotdev as σ_x = (τ/√3) · σ_y,MTOT — SigmaTau does
# not expose a ttotdev API directly.

using SigmaTau
using Printf

const HERE       = @__DIR__
const PHASE_PATH = joinpath(HERE, "..", "s32_5_12_phase.DAT")
const OUT_CSV    = joinpath(HERE, "sigmatau_deviations.csv")

const AF_SHORT = [1, 2, 4, 10, 20, 40, 100, 200, 400, 1000, 2000, 4000]
const AF_LONG  = [1, 2, 4, 10, 20, 40, 100, 200, 400, 1000, 2000, 4000, 10000]

println("Loading phase record: $PHASE_PATH")
pd = read_phase(PHASE_PATH; time_col=0, value_col=1, header=10, tau0=1.0)
println("  N = $(length(pd.x)) samples, tau0 = $(pd.tau0)s")

println("\nComputing deviations (calc_ci=true, confidence=0.683)...")

t0 = time()
results = Dict{String,Any}()
results["adev"]    = adev(pd,    AF_SHORT; calc_ci=true)
results["mdev"]    = mdev(pd,    AF_SHORT; calc_ci=true)
results["tdev"]    = tdev(pd,    AF_SHORT; calc_ci=true)
results["hdev"]    = hdev(pd,    AF_SHORT; calc_ci=true)
results["totdev"]  = totdev(pd,  AF_LONG;  calc_ci=true, correct_bias=true)
results["mtotdev"] = mtotdev(pd, AF_LONG;  calc_ci=true, correct_bias=false)
results["htotdev"] = htotdev(pd, AF_LONG;  calc_ci=true, correct_bias=true)
println("  core deviations done in $(round(time()-t0, digits=2))s")

# Derive TTOTDEV from MTOTDEV: σ_x,TTOT = (τ/√3) · σ_y,MTOT
let r = results["mtotdev"]
    f = r.tau ./ sqrt(3.0)
    results["ttotdev"] = StabilityResult(
        :ttotdev, r.tau, r.dev .* f, r.noise_type,
        r.ci_lower .* f, r.ci_upper .* f, r.edf,
    )
end

println("\nWriting CSV: $OUT_CSV")
open(OUT_CSV, "w") do io
    println(io, "deviation,AF,tau,dev,ci_lower,ci_upper,alpha,edf")
    af_for = name -> (name in ("totdev", "mtotdev", "ttotdev", "htotdev")) ? AF_LONG : AF_SHORT
    # canonical print order
    for name in ("adev","mdev","tdev","hdev","totdev","mtotdev","ttotdev","htotdev")
        r = results[name]
        af = af_for(name)
        # noise_type symbols from identify_noise → Stable32 alpha integer
        alpha_int(sym) = sym === :WHPM ?  2 :
                          sym === :FLPM ?  1 :
                          sym === :WHFM ?  0 :
                          sym === :FLFM ? -1 :
                          sym === :RWFM ? -2 : 99
        for (i, m) in enumerate(af)
            a = isempty(r.noise_type) ? 99 : alpha_int(r.noise_type[i])
            edf_val = isempty(r.edf) ? NaN : r.edf[i]
            @printf(io, "%s,%d,%.6e,%.10e,%.10e,%.10e,%d,%.6f\n",
                    name, m, r.tau[i], r.dev[i],
                    r.ci_lower[i], r.ci_upper[i], a, edf_val)
        end
    end
end
println("Done.")
