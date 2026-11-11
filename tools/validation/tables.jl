# tables.jl — rewrites the tables at each averaging time on
# docs/src/explanation/validation.md from the 8192-point record and the
# Stable32 and allantools tables stored with it.
#
#   julia --project=. tools/validation/tables.jl
#
# Each table sits under a `### <name>` heading, named as in the Stable32
# file. The text around the tables is left as it is.
using SigmaTau, Printf

const FIX  = joinpath(@__DIR__, "..", "..", "test", "fixtures", "validation")
const PAGE = joinpath(@__DIR__, "..", "..", "docs", "src", "explanation", "validation.md")

# name => (estimator, keywords, has a bias correction)
const TABLES = [
    "Overlapping Allan"    => (adev,    (;),                    false),
    "Modified Allan"       => (mdev,    (;),                    false),
    "Time"                 => (tdev,    (;),                    false),
    "Overlapping Hadamard" => (hdev,    (; remove_drift=false), false),
    "Total"                => (totdev,  (;),                    true),
    "Modified Total"       => (mtotdev, (;),                    true),
    "Hadamard Total"       => (htotdev, (; remove_drift=false), true),
]

# Rows of one estimator: averaging factor => the remaining columns.
function reference(file, name)
    rows = Dict{Int,Vector{Float64}}()
    for ln in Iterators.drop(eachline(joinpath(FIX, file)), 1)
        f = split(ln, ',')
        f[1] == name && (rows[parse(Int, f[2])] = parse.(Float64, f[3:end]))
    end
    return rows
end

s4(v) = @sprintf("%.4e", v)
s5(v) = @sprintf("%.5e", v)
band(lo, hi) = "[$(s4(lo)), $(s4(hi))]"

function table(p, name, (est, kw, biased))
    s32 = reference("stable32out/stable32_data_full.csv", name)       # τ, N, α, min, σ, max
    at  = reference("allantools_out/allantools_data_full.csv", name)  # τ, N, σ, err
    ms  = sort!(collect(keys(s32)))
    r   = est(p, ms; kw...)
    raw = biased ? est(p, ms; kw..., correct_bias=false) : r
    # Stable32 prints MTOTDEV uncorrected, so its limits go with the uncorrected call.
    shown = est === mtotdev ? raw : r
    io = IOBuffer()
    if biased
        println(io, "| τ (s) | Stable32 | allantools | Uncorrected | Default | Stable32 interval | SigmaTau interval |")
        println(io, "|---:|:---|:---|:---|:---|:---|:---|")
    else
        println(io, "| τ (s) | Stable32 | allantools | SigmaTau | Stable32 interval | SigmaTau interval |")
        println(io, "|---:|:---|:---|:---|:---|:---|")
    end
    for (k, m) in enumerate(ms)
        lo, sigma, hi = s32[m][4:6]
        mine = biased ? "$(s5(raw.dev[k])) | $(s5(r.dev[k]))" : s5(r.dev[k])
        flag = SigmaTau._alpha_from_noise(r.noise_type[k]) == s32[m][3] ? "" : "¹"   # noise type differs
        println(io, "| $m$flag | $(s4(sigma)) | $(s5(at[m][3])) | $mine | ",
                    band(lo, hi), " | ", band(shown.ci[k]...), " |")
    end
    return String(take!(io))
end

function main()
    p = read_phase(joinpath(FIX, "stable32gen.DAT"))
    page = read(PAGE, String)
    for (name, spec) in TABLES
        block = Regex("(^### \\Q$name\\E\\n(?:(?!\\|).*\\n)*?)(?:\\|.*\\n)+", "m")
        occursin(block, page) || error("no table under \"### $name\" in $PAGE")
        page = replace(page, block => SubstitutionString("\\1" * table(p, name, spec)))
    end
    write(PAGE, page)
end

main()
