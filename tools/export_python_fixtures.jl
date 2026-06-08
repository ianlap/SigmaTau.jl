# tools/export_python_fixtures.jl — emit golden parity fixtures for the Python port.
#
# Runs the MVP deviations (adev/mdev/tdev/hdev/mhdev/htdev) at ci=false on a set
# of fixed input records and writes, into the sibling `sigmatau` repo:
#   - tests/fixtures/julia/inputs/<name>.txt   single-column raw samples (so the
#       Python side reads the IDENTICAL input — RNGs differ across languages, so
#       inputs are shared, never regenerated)
#   - tests/fixtures/julia/julia_reference.csv long-format reference outputs
#       (input, data_kind, deviation, grid, m, tau, dev)
#
# Run from the SigmaTau.jl repo root:
#     julia --project=. tools/export_python_fixtures.jl [path/to/sigmatau]

using SigmaTau
using Printf
using Random

const OUT_ROOT = length(ARGS) >= 1 ? ARGS[1] :
                 normpath(joinpath(@__DIR__, "..", "..", "sigmatau"))
const FIX_DIR  = joinpath(OUT_ROOT, "tests", "fixtures", "julia")
const IN_DIR   = joinpath(FIX_DIR, "inputs")
mkpath(IN_DIR)

const DEVS = (
    ("adev", adev), ("mdev", mdev), ("tdev", tdev),
    ("hdev", hdev), ("mhdev", mhdev), ("htdev", htdev),
    ("totdev", totdev), ("mtie", mtie), ("pdev", pdev),
)

# The modified-total family is loop-heavy in pure Python, so it is exported only
# on the small synthetic records (octave grid) — enough to validate parity
# without making the Python test suite slow.
const MODTOTAL = (
    ("mtotdev", mtotdev), ("ttotdev", ttotdev),
    ("htotdev", htotdev), ("mhtotdev", mhtotdev),
)

# Emit one CSV row per τ, with the full ci=true output (noise type, EDF, χ² CI
# bounds). Deviations without a CI model (mtie) leave those columns empty/NaN.
function write_rows(io, name, kind, dname, gname, m, r)
    for i in eachindex(m)
        nt  = isempty(r.noise_type) ? "" : string(r.noise_type[i])
        edf = isempty(r.edf)        ? NaN : r.edf[i]
        lo  = isempty(r.ci_lower)   ? NaN : r.ci_lower[i]
        hi  = isempty(r.ci_upper)   ? NaN : r.ci_upper[i]
        @printf(io, "%s,%s,%s,%s,%d,%.17e,%.17e,%s,%.17e,%.17e,%.17e\n",
                name, kind, dname, gname, m[i], r.tau[i], r.dev[i], nt, edf, lo, hi)
    end
end

# Write a single-column sample file (full Float64 precision).
function write_samples(path, v)
    open(path, "w") do io
        for s in v
            @printf(io, "%.17e\n", s)
        end
    end
end

# Fixed input records. Content is arbitrary for parity (Python reads the same
# file); we vary character to exercise the kernels. tau0 = 1.0 throughout.
function build_records()
    recs = Tuple{String,Symbol,Any}[]

    # 1. The real Stable32 composite phase record (10-line header, 8192 samples).
    s32 = joinpath(@__DIR__, "..", "reference", "validation", "stable32gen.DAT")
    if isfile(s32)
        lines = readlines(s32)
        x = parse.(Float64, strip.(lines[11:end]))
        write_samples(joinpath(IN_DIR, "stable32gen_phase.txt"), x)
        push!(recs, ("stable32gen_phase", :phase, PhaseData(x, 1.0)))
    else
        @warn "stable32gen.DAT not found; skipping that record" path=s32
    end

    # 2. Synthetic random-walk phase (N=1024).
    rng = MersenneTwister(20260607)
    xs = cumsum(randn(rng, 1024)) .* 1e-9
    write_samples(joinpath(IN_DIR, "synth_phase.txt"), xs)
    push!(recs, ("synth_phase", :phase, PhaseData(xs, 1.0)))

    # 3. Synthetic white frequency (N=1024).
    ys = randn(rng, 1024) .* 1e-12
    write_samples(joinpath(IN_DIR, "synth_freq.txt"), ys)
    push!(recs, ("synth_freq", :frequency, FrequencyData(ys, 1.0)))

    return recs
end

function main()
    recs = build_records()
    csv = joinpath(FIX_DIR, "julia_reference.csv")
    open(csv, "w") do io
        println(io, "input,data_kind,deviation,grid,m,tau,dev,noise_type,edf,ci_lower,ci_upper")
        for (name, kind, data) in recs
            n = kind === :phase ? length(data.x) : length(data.y)
            # Octave for every record; add the dense all-tau grid for the small
            # synthetic ones (the big record's all-tau grid is needlessly large).
            grids = n <= 2048 ? (("octave", Octave), ("alltaus", AllTaus)) :
                                (("octave", Octave),)
            for (gname, gmode) in grids
                for (dname, fn) in DEVS
                    m = tau_values(gmode, n, Symbol(dname))
                    r = fn(data, m)  # ci=true defaults (oracle behavior)
                    write_rows(io, name, kind, dname, gname, m, r)
                end

                # Modified-total family: synthetic records, octave grid only.
                if gname == "octave" && n <= 2048
                    for (dname, fn) in MODTOTAL
                        m = tau_values(gmode, n, Symbol(dname))
                        r = fn(data, m)  # ci=true, correct_bias=true defaults
                        write_rows(io, name, kind, dname, gname, m, r)
                    end
                end
            end
        end
    end

    # Spectral estimators (Sy / Sx / L) — separate CSV (freq/psd shaped, not
    # deviation rows). Default Welch params; L at a fixed 10 MHz carrier.
    const_fc = 1.0e7
    scsv = joinpath(FIX_DIR, "spectral_reference.csv")
    open(scsv, "w") do io
        println(io, "input,estimator,f_carrier,idx,freq,psd")
        for (name, kind, data) in recs
            kind === :phase || kind === :frequency || continue
            (kind === :phase && length(data.x) > 2048) && continue  # synth only
            (kind === :frequency && length(data.y) > 2048) && continue
            for (ename, r) in (("Sy", Sy(data)), ("Sx", Sx(data)),
                               ("L", L(data; f_carrier=const_fc)))
                fc = ename == "L" ? const_fc : NaN
                for i in eachindex(r.freq)
                    @printf(io, "%s,%s,%.6e,%d,%.17e,%.17e\n",
                            name, ename, fc, i - 1, r.freq[i], r.psd[i])
                end
            end
        end
    end

    # IO fixtures: a fillgaps parity case + Julia-saved result/suite TSVs (for
    # the Python port's cross-language load compatibility).
    synth = first(d for (nm, _, d) in recs if nm == "synth_phase")

    xgap = copy(synth.x)
    for r in (101:106, 301:301, 700:710, 1000:1002)
        xgap[r] .= NaN
    end
    filled = fillgaps(PhaseData(xgap, 1.0)).x
    open(joinpath(FIX_DIR, "fillgaps_reference.csv"), "w") do io
        println(io, "idx,input,filled")
        for i in eachindex(xgap)
            instr = isnan(xgap[i]) ? "NaN" : @sprintf("%.17e", xgap[i])
            @printf(io, "%d,%s,%.17e\n", i - 1, instr, filled[i])
        end
    end

    save_result(joinpath(FIX_DIR, "adev_result.tsv"), adev(synth))
    save_suite(joinpath(FIX_DIR, "suite.tsv"), stability(synth))

    println("Wrote fixtures to $FIX_DIR")
end

main()
