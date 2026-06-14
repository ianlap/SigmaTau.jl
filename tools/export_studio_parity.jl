#!/usr/bin/env julia
# Regenerate SigmaTauStudio's MHTOTDEV parity fixtures from the CURRENT SigmaTau.jl oracle.
#
#   julia --project=. tools/export_studio_parity.jl
#
# Writes into ../../SigmaTauStudio/crates/sigmatau/tests/fixtures (sibling repo):
#   * julia_reference.csv  — the `mhtotdev,*` rows are replaced with current-oracle
#                            values (every other estimator's rows are left untouched).
#   * mhtot_<name>.dat      — 6 seeded N=8192 phase records: one per noise type
#                            (whpm/flpm/whfm/flfm/rwfm → α=2,1,0,−1,−2) + one `drifted`.
#   * mhtotdev_parity.csv   — two layers per record: `raw` (correct_bias=false,
#                            remove_drift=false) and `full` (the default path).
#
# The dedicated fixture exists because a single record cannot exercise all five
# per-α bias/EDF rows, and never tests the wrapper's drift removal. Inputs are
# committed verbatim so the Rust engine reads identical samples (parity does not
# depend on cross-language RNG matching).

using SigmaTau
using Random
using Printf

const FIX = abspath(joinpath(@__DIR__, "..", "..", "SigmaTauStudio",
                             "crates", "sigmatau", "tests", "fixtures"))

const ALPHA = Dict(:WHPM => 2, :FLPM => 1, :WHFM => 0, :FLFM => -1, :RWFM => -2)
alpha_of(sym::Symbol) = get(ALPHA, sym, 0)

g(x::Real) = @sprintf("%.17e", x)   # match the existing CSV's full-precision format

# ── Part A: replace the mhtotdev rows in julia_reference.csv ────────────────
function read_record_dat(path)
    x = Float64[]
    indata = false
    for line in eachline(path)
        s = strip(line)
        if startswith(s, "# Header End"); indata = true; continue; end
        indata && !isempty(s) && push!(x, parse(Float64, s))
    end
    return x
end

function mhtot_reference_rows(pd)
    r = mhtotdev(pd; ci=true, correct_bias=true, remove_drift=true)
    rows = String[]
    for i in eachindex(r.tau)
        m = round(Int, r.tau[i] / pd.tau0)
        push!(rows, join(["mhtotdev", string(m), g(r.tau[i]), g(r.dev[i]),
                          g(r.ci[i].lo), g(r.ci[i].hi), string(alpha_of(r.noise_type[i])),
                          g(r.edf[i]), string(r.neff[i])], ","))
    end
    return rows
end

function regen_reference!()
    path = joinpath(FIX, "julia_reference.csv")
    lines = readlines(path)
    header, body = lines[1], lines[2:end]
    kept = filter(l -> !startswith(l, "mhtotdev,"), body)
    pd = PhaseData(read_record_dat(joinpath(FIX, "record.dat")), 1.0)
    newrows = mhtot_reference_rows(pd)
    open(path, "w") do io
        println(io, header)
        foreach(l -> println(io, l), kept)
        foreach(l -> println(io, l), newrows)
    end
    println("Part A: replaced $(length(newrows)) mhtotdev rows in julia_reference.csv")
    println("        m=1 σ=$(g(parse(Float64, split(newrows[1], ',')[4])))")
end

# ── Part B: dedicated MHTOTDEV fixture (5 noise types + 1 drifted) ──────────
function gen_records()
    recs = Tuple{String,Vector{Float64}}[]
    seed0 = 20260614
    for (i, nt) in enumerate((:WHPM, :FLPM, :WHFM, :FLFM, :RWFM))
        pd = noise_gen(PhaseData, 8192, 1.0;
                       sigma1=Dict(ALPHA[nt] => 1e-12), rng=Xoshiro(seed0 + i))
        push!(recs, (lowercase(string(nt)), pd.x))
    end
    # drifted: clean WHFM + a deterministic frequency drift D (linear in freq ⇒
    # quadratic in phase), large enough to dominate, so remove_drift is decisive.
    clean = noise_gen(PhaseData, 8192, 1.0; sigma1=Dict(0 => 1e-12), rng=Xoshiro(seed0 + 99))
    D = 1e-13
    drifted = [clean.x[i] + 0.5 * D * ((i - 1) * 1.0)^2 for i in eachindex(clean.x)]
    push!(recs, ("drifted", drifted))
    return recs
end

function write_layer(io, name, layer, r; full::Bool)
    for i in eachindex(r.tau)
        m = round(Int, r.tau[i])
        if full
            println(io, join([name, layer, string(m), g(r.tau[i]), g(r.dev[i]),
                              g(r.ci[i].lo), g(r.ci[i].hi), string(alpha_of(r.noise_type[i])),
                              g(r.edf[i]), string(r.neff[i])], ","))
        else
            println(io, join([name, layer, string(m), g(r.tau[i]), g(r.dev[i]),
                              "nan", "nan", "0", "nan", string(r.neff[i])], ","))
        end
    end
end

function write_fixture()
    recs = gen_records()
    for (name, x) in recs
        open(joinpath(FIX, "mhtot_$(name).dat"), "w") do io
            foreach(v -> println(io, g(v)), x)
        end
    end
    open(joinpath(FIX, "mhtotdev_parity.csv"), "w") do io
        println(io, "record,layer,m,tau,sigma,ci_lo,ci_hi,alpha,edf,neff")
        for (name, x) in recs
            pd = PhaseData(copy(x), 1.0)
            write_layer(io, name, "raw",
                        mhtotdev(pd; ci=false, correct_bias=false, remove_drift=false); full=false)
            write_layer(io, name, "full",
                        mhtotdev(pd; ci=true, correct_bias=true, remove_drift=true); full=true)
        end
    end
    println("Part B: wrote mhtotdev_parity.csv + $(length(recs)) record files")
end

regen_reference!()
write_fixture()
println("done → $FIX")
