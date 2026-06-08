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
)

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
        println(io, "input,data_kind,deviation,grid,m,tau,dev")
        for (name, kind, data) in recs
            n = kind === :phase ? length(data.x) : length(data.y)
            # Octave for every record; add the dense all-tau grid for the small
            # synthetic ones (the big record's all-tau grid is needlessly large).
            grids = n <= 2048 ? (("octave", Octave), ("alltaus", AllTaus)) :
                                (("octave", Octave),)
            for (gname, gmode) in grids
                for (dname, fn) in DEVS
                    m = tau_values(gmode, n, Symbol(dname))
                    r = fn(data, m; ci=false)
                    for i in eachindex(m)
                        @printf(io, "%s,%s,%s,%s,%d,%.17e,%.17e\n",
                                name, kind, dname, gname, m[i], r.tau[i], r.dev[i])
                    end
                end
            end
        end
    end
    println("Wrote fixtures to $FIX_DIR")
end

main()
