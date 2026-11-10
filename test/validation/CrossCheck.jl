# CrossCheck.jl — loaders and comparison drivers for the reference cross-check
# suite. Used by test/validation/runtests.jl (assertions) and by
# tools/validation/report.jl (the Markdown discrepancy table), so both read the
# same rows.
#
# Every comparison produces `Row`s: one reference value, SigmaTau's value, the
# relative error, the tolerance it is held to and a short note naming the
# convention when the two programs report different statistics on purpose.
module CrossCheck

using SigmaTau
using SigmaTau: calculate_edf, bias_correction, confidence_intervals, identify_noise,
                _noise_id_lag1acf, _calc_edf_core, _freq_to_phase

export Row, FIX, ALPHAS, noise_record, phasedat, stable32gen, nbs14_frequency,
       all_rows, relerr

const FIX = normpath(joinpath(@__DIR__, "..", "fixtures", "validation"))
const AT  = get(ENV, "SIGMATAU_ALLANTOOLS_DIR", joinpath(FIX, "allantools"))
const S32 = joinpath(FIX, "stable32_phasedat")
const ALPHAS = (2, 1, 0, -1, -2, -3, -4)

struct Row
    group::String       # e.g. "allantools σ", "Stable32 CI"
    dataset::String
    estimator::String
    m::Int
    ref::Float64
    got::Float64
    tol::Float64        # relative tolerance this row is held to
    note::String
    floor::Float64      # error is relative to max(|ref|, floor)
end
Row(g, d, e, m, ref, got, tol, note) = Row(g, d, e, m, ref, got, tol, note, floatmin())

relerr(r::Row) = abs(r.got - r.ref) / max(abs(r.ref), r.floor)
# tol = Inf marks a row that is reported but not asserted (a documented
# convention or a reference that has no value there).
passes(r::Row) = isinf(r.tol) || (isnan(r.ref) && isnan(r.got)) || relerr(r) <= r.tol

# ── data ────────────────────────────────────────────────────────────────────

_floats(path; header_end=false) = begin
    out = Float64[]
    started = !header_end
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        if !started
            startswith(s, "# Header End") && (started = true)
            continue
        end
        startswith(s, "#") && continue
        push!(out, parse(Float64, last(split(s))))
    end
    out
end

alpha_file(a) = a > 0 ? "alpha_p$a.txt" : a == 0 ? "alpha_0.txt" : "alpha_m$(-a).txt"
noise_record(a) = PhaseData(_floats(joinpath(FIX, "noise", alpha_file(a))), 1.0)
phasedat()      = PhaseData(_floats(joinpath(S32, "PHASE.DAT")), 1.0)
stable32gen()   = PhaseData(_floats(joinpath(FIX, "stable32gen.DAT"); header_end=true), 1.0)

"NBS 1000-point frequency set (SP1065 §12.3): Lehmer LCG, seed 1234567890."
function nbs14_frequency()
    n = Vector{Int}(undef, 1000)
    n[1] = 1234567890
    for i in 1:999
        n[i+1] = mod(16807 * n[i], 2147483647)
    end
    return n ./ 2147483647
end

function dataset(name::AbstractString)
    name == "phasedat"    && return phasedat()
    name == "stable32gen" && return stable32gen()
    startswith(name, "alpha_") && return noise_record(parse(Int, name[7:end]))
    error("unknown dataset $name")
end

function _csv(path)
    lines = readlines(path)
    hdr = split(lines[1], ',')
    rows = [Dict(zip(hdr, split(l, ','))) for l in lines[2:end] if !isempty(l)]
    return rows
end
_f(s) = s == "nan" ? NaN : parse(Float64, s)

# ── allantools σ ────────────────────────────────────────────────────────────

# SigmaTau call reproducing each allantools estimator. allantools reports the
# raw (bias-free) total estimators, so those run with correct_bias=false.
const AT_MAP = Dict(
    "oadev"   => (p, m) -> adev(p, [m]; ci=false).dev[1],
    "mdev"    => (p, m) -> mdev(p, [m]; ci=false).dev[1],
    "tdev"    => (p, m) -> tdev(p, [m]; ci=false).dev[1],
    "ohdev"   => (p, m) -> hdev(p, [m]; ci=false).dev[1],
    "totdev"  => (p, m) -> totdev(p, [m]; ci=false, correct_bias=false).dev[1],
    "mtotdev" => (p, m) -> mtotdev(p, [m]; ci=false, correct_bias=false).dev[1],
    "ttotdev" => (p, m) -> ttotdev(p, [m]; ci=false, correct_bias=false).dev[1],
    "htotdev" => (p, m) -> htotdev(p, [m]; ci=false, correct_bias=false).dev[1],
    "tierms"  => (p, m) -> tierms(p, [m]).dev[1],
    "mtie"    => (p, m) -> mtie(p, [m]).dev[1],
    "pdev"    => (p, m) -> pdev(p, [m]; ci=false).dev[1],
    "theo1"   => (p, m) -> theo1(p, [m]; ci=false, correct_bias=false).dev[1],
)

# Per-estimator tolerance (relative) and the note that explains it.
function at_tol(est, ds, m)
    # MDEV/TDEV run on prefix sums: on the random-walk-of-random-walk records
    # (α ≤ −3) that costs ~1e-9 relative to allantools' direct sums.
    est in ("mdev", "tdev") && return (1e-8, "")
    return (1e-9, "")
end

function allantools_dev_rows()
    rows = Row[]
    cache = Dict{String,PhaseData}()
    for r in _csv(joinpath(AT, "devs.csv"))
        est = r["estimator"]
        haskey(AT_MAP, est) || continue          # adev / hdev: non-overlapped, not in SigmaTau
        ds = r["dataset"]
        p = get!(() -> dataset(ds), cache, ds)
        m = parse(Int, r["m"])
        ref = _f(r["dev"])
        got = try AT_MAP[est](p, m) catch; NaN end
        tol, note = at_tol(est, ds, m)
        push!(rows, Row("allantools σ", ds, est, m, ref, got, tol, note))
    end
    return rows
end

# ── allantools EDF / CI / noise ID / spectra / dynamic ─────────────────────

_alpha_sym(a) = a == 2 ? :WHPM : a == 1 ? :FLPM : a == 0 ? :WHFM : a == -1 ? :FLFM : :RWFM

function allantools_edf_rows()
    rows = Row[]
    for r in _csv(joinpath(AT, "edf.csv"))
        fam = r["family"]; a = parse(Int, r["alpha"]); d = parse(Int, r["d"])
        m = parse(Int, r["m"]); N = parse(Int, r["N"]); ref = _f(r["edf"])
        if fam == "overlapping"
            got = _calc_edf_core(a, d, m, m, m, N)
            est = d == 2 ? "oadev" : "ohdev"
        elseif fam == "modified"
            got = _calc_edf_core(a, d, m, 1, m, N)
            est = d == 2 ? "mdev" : "mhdev"
        elseif fam in ("totdev", "mtotdev")
            # SP1065 Tables 7/8: edf = b·T/τ − c. allantools takes T/τ = N/m;
            # SigmaTau takes the record span T = (N − 1)·τ0 (a 1/N relative
            # difference). Compare the table formula at allantools' T here.
            est = fam
            got = calculate_edf(Symbol(fam), [1.0], [_alpha_sym(a)], [m], [Float64(m)], N, Float64(N))[1]
        else
            continue   # nonoverlapping, simple: no SigmaTau counterpart
        end
        tol, note = 1e-9, ""
        if est == "mtotdev" && m < 16
            tol, note = Inf, "SigmaTau uses MDEV's EDF below τ = 16τ0, the edge of SP1065 Table 8"
        elseif isnan(ref) && est == "ohdev" && a == -4
            tol, note = Inf, "allantools' edf_greenhall cannot reach Table 2's α = −4 row"
        elseif isnan(ref) && est == "oadev" && a == 2
            tol, note = Inf, "allantools raises for white PM when ⌈M/S⌉ ≤ d; SigmaTau sums the series"
        end
        push!(rows, Row("allantools EDF", "N=$N α=$a", est, m, ref, got, tol, note))
    end
    return rows
end

function allantools_ci_rows()
    rows = Row[]
    for r in _csv(joinpath(AT, "ci.csv"))
        e = _f(r["edf"]); c = _f(r["confidence"])
        lo, hi = confidence_intervals([1.0], [e], c)
        push!(rows, Row("allantools CI", "edf=$(r["edf"][1:6]) c=$(r["confidence"][1:6])", "lo", 0, _f(r["lo_factor"]), lo[1], 1e-9, ""))
        push!(rows, Row("allantools CI", "edf=$(r["edf"][1:6]) c=$(r["confidence"][1:6])", "hi", 0, _f(r["hi_factor"]), hi[1], 1e-9, ""))
    end
    return rows
end

function allantools_noiseid_rows()
    rows = Row[]
    cache = Dict{String,PhaseData}()
    for r in _csv(joinpath(AT, "noiseid.csv"))
        ds = r["dataset"]; p = get!(() -> dataset(ds), cache, ds)
        m = parse(Int, r["m"]); dmax = parse(Int, r["dmax"])
        alpha, a_int, d, rho = _noise_id_lag1acf(Vector{Float64}(p.x), m, 0, dmax; detrend=true)
        push!(rows, Row("allantools noise ID", ds, "α (dmax=$dmax)", m, _f(r["alpha"]), alpha, 1e-6, ""))
    end
    return rows
end

# ── Stable32 (PHASE.DAT runs shipped with allantools' test suite) ─────────

"Parse a Stable32 results file → rows of (m, n, alpha, lo, sigma, hi) or (m, n, value)."
function read_s32(file)
    rows = NamedTuple[]
    for ln in eachline(joinpath(S32, file))
        s = strip(ln)
        (isempty(s) || startswith(s, "#")) && continue
        v = parse.(Float64, split(s))
        if length(v) == 7
            push!(rows, (m=Int(v[1]), n=v[3], alpha=Int(v[4]), lo=v[5], sigma=v[6], hi=v[7]))
        elseif length(v) == 4
            push!(rows, (m=Int(v[1]), n=v[3], alpha=0, lo=0.0, sigma=v[4], hi=0.0))
        end
    end
    return rows
end

# file => (estimator label, SigmaTau call, EDF method, σ_x scale (m -> factor), note)
const S32_RUNS = [
    ("phase_dat_oadev_octave.txt",   "oadev",   (p, m) -> adev(p, m),   :adev,    m -> 1.0, ""),
    ("phase_dat_mdev_octave.txt",    "mdev",    (p, m) -> mdev(p, m),   :mdev,    m -> 1.0, ""),
    ("phase_dat_tdev_octave.txt",    "tdev",    (p, m) -> tdev(p, m),   :mdev,    m -> m / sqrt(3), ""),
    ("phase_dat_ohdev_octave.txt",   "ohdev",   (p, m) -> hdev(p, m; remove_drift=false), :hdev, m -> 1.0,
        "Stable32 identifies noise on the record as given; remove_drift=false"),
    ("phase_dat_totdev_octave.txt",  "totdev",  (p, m) -> totdev(p, m), :totdev,  m -> 1.0, ""),
    ("phase_dat_mtotdev_octave.txt", "mtotdev", (p, m) -> mtotdev(p, m; correct_bias=false), :mtotdev, m -> 1.0,
        "Stable32 reports MTOT raw; correct_bias=false"),
    ("phase_dat_ttotdev_octave.txt", "ttotdev", (p, m) -> ttotdev(p, m; correct_bias=false), :mtotdev, m -> m / sqrt(3),
        "Stable32 reports TTOT raw; correct_bias=false"),
    ("phase_dat_htotdev_octave.txt", "htotdev", (p, m) -> htotdev(p, m; remove_drift=false), :htotdev, m -> 1.0,
        "Stable32 identifies noise on the record as given; remove_drift=false"),
]

# σ-only runs over every τ (no CI columns in these files).
const S32_ALLTAU = [
    ("phase_dat_oadev.txt",  "oadev",  (p, m) -> adev(p, m; ci=false)),
    ("phase_dat_mdev.txt",   "mdev",   (p, m) -> mdev(p, m; ci=false)),
    ("phase_dat_tdev.txt",   "tdev",   (p, m) -> tdev(p, m; ci=false)),
    ("phase_dat_ohdev.txt",  "ohdev",  (p, m) -> hdev(p, m; ci=false)),
    ("phase_dat_totdev.txt", "totdev", :totdev_s32alpha),
    ("phase_dat_tierms.txt", "tierms", (p, m) -> tierms(p, m)),
    ("phase_dat_mtie.txt",   "mtie",   (p, m) -> mtie(p, m)),
    ("phase_dat_theo1_alpha0.txt",        "theo1", (p, m) -> theo1(p, m; ci=false, correct_bias=false)),
    ("phase_dat_theo1_alpha0_decade.txt", "theo1", (p, m) -> theo1(p, m; ci=false, correct_bias=false)),
]

const S32_TOL_SIGMA = 1e-4   # Stable32 prints 5 significant figures
const S32_TOL_CI    = 1e-3

const NOISE_ID_NOTE = "noise ID differs: N/m < 30 is the B1/R(n) fallback, which reads 7 averaged samples here"

# CI rows outside the published EDF models: reported, not asserted.
function s32_ci_tol(est, m)
    if est in ("mtotdev", "ttotdev") && m >= 16
        return (Inf, "Stable32's MTOT EDF at τ = 16τ0 is not SP1065 Table 8 (implied EDF 81 vs 68)")
    elseif est == "htotdev" && m < 16
        return (Inf, "below τ = 16τ0 the FCS 2001 HTOT EDF table does not apply; Stable32's fallback is undocumented")
    end
    return (S32_TOL_CI, "")
end

function stable32_rows()
    p = phasedat()
    N = length(p.x); T = (N - 1) * p.tau0
    rows = Row[]
    for (file, est, f, method, scale, note) in S32_RUNS
        s = read_s32(file)
        ms = [r.m for r in s]
        res = f(p, ms)
        for (k, r) in enumerate(s)
            push!(rows, Row("Stable32 σ", "PHASE.DAT", est, r.m, r.sigma, res.dev[k], S32_TOL_SIGMA, note))
            a_st = SigmaTau._alpha_from_noise(res.noise_type[k])
            same = a_st == r.alpha
            idnote = same ? "" : NOISE_ID_NOTE
            push!(rows, Row("Stable32 α", "PHASE.DAT", est, r.m, r.alpha, a_st, same ? 0.0 : Inf, idnote))
            # CI with SigmaTau's own noise ID (the default API output)
            tol, cinote = s32_ci_tol(est, r.m)
            same || ((tol, cinote) = (Inf, idnote))
            push!(rows, Row("Stable32 CI", "PHASE.DAT", est * " lo", r.m, r.lo, res.ci[k].lo, tol, cinote))
            push!(rows, Row("Stable32 CI", "PHASE.DAT", est * " hi", r.m, r.hi, res.ci[k].hi, tol, cinote))
        end
    end
    for (file, est, f) in S32_ALLTAU
        s = read_s32(file)
        ms = [r.m for r in s]
        if f === :totdev_s32alpha
            # Stable32 identifies α = 0 or +1 out to m = 1000 here (N/m < 30 is the
            # B1 fallback, where SigmaTau finds FLFM/RWFM at m ≳ 210 and applies
            # the TOTVAR bias). Apply the bias for Stable32's α so the rows
            # compare the estimator, not the noise ID.
            raw = totdev(p, ms; ci=false, correct_bias=false).dev
            B = bias_correction([_alpha_sym(r.alpha) for r in s], :totvar, Float64.(ms), T)
            devs = raw ./ sqrt.(B)
            note = "bias for Stable32's α (noise ID differs at N/m < 5)"
        else
            devs = f(p, ms).dev
            note = ""
        end
        for (k, r) in enumerate(s)
            push!(rows, Row("Stable32 σ (all τ)", "PHASE.DAT", est, r.m, r.sigma, devs[k], S32_TOL_SIGMA, note))
        end
    end
    return rows
end

"Stable32 CI recomputed with Stable32's own α, isolating the EDF model from noise ID."
function stable32_ci_forced_alpha_rows()
    p = phasedat()
    N = length(p.x); T = (N - 1) * p.tau0
    rows = Row[]
    for (file, est, f, method, scale, note) in S32_RUNS
        for r in read_s32(file)
            m = r.m
            σ = r.sigma / scale(m)                       # back to the σ_y statistic
            edf = calculate_edf(method, [σ], [_alpha_sym(r.alpha)], [m], [Float64(m)], N, T)[1]
            lo, hi = confidence_intervals([σ], [edf], 0.683)
            tol, cinote = s32_ci_tol(est, m)
            push!(rows, Row("Stable32 CI (S32 α)", "PHASE.DAT", est * " lo", m, r.lo, lo[1] * scale(m), tol, cinote))
            push!(rows, Row("Stable32 CI (S32 α)", "PHASE.DAT", est * " hi", m, r.hi, hi[1] * scale(m), tol, cinote))
        end
    end
    return rows
end

# ── NBS 1000-point published values (SP1065 §12.3; allantools nbs14 test) ──

# (estimator, call, m values, values, note). m = τ/τ0 on the frequency record.
const NBS14_TABLE = [
    ("oadev",   (d, m) -> adev(d, m; ci=false),   [1, 10, 100], [2.922319e-01, 9.159953e-02, 3.241343e-02], ""),
    ("mdev",    (d, m) -> mdev(d, m; ci=false),   [1, 10, 100], [2.922319e-01, 6.172376e-02, 2.170921e-02], ""),
    ("tdev",    (d, m) -> tdev(d, m; ci=false),   [1, 10, 100], [1.687202e-01, 3.563623e-01, 1.253382e+00], ""),
    ("ohdev",   (d, m) -> hdev(d, m; ci=false),   [1, 10, 100], [2.943883e-01, 9.581083e-02, 3.237638e-02], ""),
    ("totdev",  (d, m) -> totdev(d, m),           [1, 10, 100], [2.922319e-01, 9.134743e-02, 3.406530e-02], ""),
    ("htotdev", (d, m) -> htotdev(d, m; correct_bias=false).dev ./ [mm == 1 ? 1.0 : sqrt(0.995) for mm in m],
        [1, 10, 100], [2.943883e-01, 9.614787e-02, 3.058103e-02],
        "published values carry the white-FM TotHvar bias (B = 0.995) above m = 1; the set is white FM by construction"),
    ("mtotdev", (d, m) -> mtotdev(d, m; correct_bias=false), [1, 10, 100], [2.0664e-01, 5.5529e-02, 1.9547e-02],
        "Stable32 1.60 raw run (SP1065 Table 31's 2.418528e-01 is the Stable32 1.53 algorithm)"),
    ("ttotdev", (d, m) -> ttotdev(d, m; correct_bias=false), [1, 10, 100], [1.1930e-01, 3.2060e-01, 1.1285e+00],
        "Stable32 1.60 raw run"),
    ("theo1",   (d, m) -> theo1(d, m; ci=false, correct_bias=false), [10, 100, 1000], [1.0757e-01, 3.1789e-02, 5.0524e-03],
        "Stable32 raw Thêo1 at AF = m (τ = 0.75 m)"),
    ("pdev",    (d, m) -> pdev(d, m; ci=false),   [1, 2, 4, 8, 16, 32, 64, 128, 256],
        [2.9223187810675200e-01, 2.1445233564252639e-01, 1.5618112158618463e-01, 1.1709745745448434e-01,
         6.9029585189839343e-02, 4.9749707730398392e-02, 3.8947417330713739e-02, 3.0862392741372108e-02,
         1.2447414341332683e-02], "SigmaTheta Pdev"),
]

function nbs14_rows()
    fd = FrequencyData(nbs14_frequency(), 1.0)
    rows = Row[]
    for (est, f, ms, vals, note) in NBS14_TABLE
        res = f(fd, ms)
        got = res isa StabilityResult ? res.dev : res
        tol = est == "pdev" ? 1e-11 : (length(note) > 0 && occursin("1.60", note)) || est == "theo1" ? 1e-4 : 2e-6
        for (m, v, g) in zip(ms, vals, got)
            push!(rows, Row("NBS 1000-point", "NBS14", est, m, v, g, tol, note))
        end
    end
    return rows
end

# ── scipy Welch spectra and per-window allantools (dynamic maps) ────────────

function welch_rows()
    rows = Row[]
    cache = Dict{String,PhaseData}()
    byk = Dict{Tuple{String,String,Int},Vector{Float64}}()
    for r in _csv(joinpath(AT, "welch.csv"))
        push!(get!(byk, (r["dataset"], r["kind"], parse(Int, r["nperseg"])), Float64[]), _f(r["psd"]))
    end
    for ((ds, kind, nps), ref) in sort(collect(byk); by=first)
        p = get!(() -> dataset(ds), cache, ds)
        res = kind == "Sy" ? psd_y(p; nperseg=nps) : psd_x(p; nperseg=nps)
        # Error relative to the spectral peak: FFT round-off scales with the
        # total power, and a random-run phase record spans ~20 decades of S_x.
        pk = maximum(ref)
        for (k, (a, b)) in enumerate(zip(ref, res.psd))
            push!(rows, Row("scipy Welch", ds, kind, k - 1, a, b, 1e-9, "error relative to the peak", pk))
        end
    end
    return rows
end

function dynamic_rows()
    rows = Row[]
    cache = Dict{Tuple{String,String},DynamicStabilityResult}()
    ms = [1, 2, 4, 16, 64, 256]
    for r in _csv(joinpath(AT, "dynamic.csv"))
        ds = r["dataset"]; kind = r["kind"]
        res = get!(cache, (ds, kind)) do
            p = dataset(ds)
            kind == "dadev" ? dadev(p, ms; window=1024, step=512) : dhdev(p, ms; window=1024, step=512)
        end
        i = (parse(Int, r["start"]) - 1) ÷ 512 + 1
        j = findfirst(==(parse(Int, r["m"])), ms)
        push!(rows, Row("allantools per-window", ds, kind, ms[j], _f(r["dev"]), res.dev[i, j], 1e-9, "window $i"))
    end
    return rows
end

# ── allantools gradev on records with gaps ──────────────────────────────────

# Samples blanked before both programs run; the same spans are in
# tools/validation/gen_allantools.py.
const GAP_SPANS = ((101, 160), (777, 777), (1500, 1503), (3001, 3100))

function gap_rows()
    rows = Row[]
    cache = Dict{String,StabilityResult}()
    for r in _csv(joinpath(AT, "gaps.csv"))
        ds = r["dataset"]
        res = get!(cache, ds) do
            x = copy(dataset(ds).x)
            for (lo, hi) in GAP_SPANS
                x[lo:hi] .= NaN
            end
            ms = [parse(Int, q["m"]) for q in _csv(joinpath(AT, "gaps.csv")) if q["dataset"] == ds]
            adev(PhaseData(x, 1.0), ms; ci=false)
        end
        m = parse(Int, r["m"])
        j = findfirst(==(m), Int.(res.tau))
        push!(rows, Row("allantools gradev", ds, "adev", m, _f(r["dev"]), res.dev[j], 1e-12, ""))
        push!(rows, Row("allantools gradev", ds, "terms", m, parse(Float64, r["n"]), Float64(res.neff[j]), 0.0, ""))
    end
    return rows
end

# ── scipy least-squares drift fits ──────────────────────────────────────────

function drift_rows()
    rows = Row[]
    y0 = diff(noise_record(0).x)
    t = collect(0.0:length(y0)-1)
    model_y = Dict(
        "linear"    => (a, b, c) -> a .+ b .* t,
        "quadratic" => (a, b, c) -> a .+ b .* t .+ c .* t .^ 2,
        "log"       => (a, b, c) -> a .* log1p.(b .* t) .+ c,
        "diffusion" => (a, b, c) -> a .+ b .* sqrt.(t .+ c))
    cache = Dict{String,DriftFit}()
    for r in _csv(joinpath(AT, "drift.csv"))
        model = r["model"]
        f = get!(cache, model) do
            inj = (_f(r["inj_a"]), _f(r["inj_b"]), _f(r["inj_c"]))
            drift_fit(FrequencyData(y0 .+ model_y[model](inj...), 1.0); model=Symbol(model))
        end
        k = findfirst(r["param"], "abc")[1]
        # The nonlinear residual is flat to rounding near its minimum, so
        # their parameters agree to ~1e-8 rather than to the last digit.
        tol = model in ("linear", "quadratic") ? 1e-12 : 1e-7
        push!(rows, Row("scipy least squares", "alpha_+0", model, 0, _f(r["value"]), f.coef[k], tol, r["param"]))
        push!(rows, Row("scipy least squares", "alpha_+0", model, 0, _f(r["stderr"]), f.stderr[k], tol, "σ(" * r["param"] * ")"))
    end
    return rows
end

"Every comparison row, in report order."
all_rows() = vcat(allantools_dev_rows(), nbs14_rows(), stable32_rows(), stable32_ci_forced_alpha_rows(),
                  allantools_edf_rows(), allantools_ci_rows(), allantools_noiseid_rows(),
                  welch_rows(), dynamic_rows(), gap_rows(), drift_rows())

end # module
