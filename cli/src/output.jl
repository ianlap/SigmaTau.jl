# JSON and CSV, written by hand (no JSON package: it would not trim).

function json_string(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"';      print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\r'; print(io, "\\r")
        elseif c == '\t'; print(io, "\\t")
        elseif c < ' ';   print(io, "\\u00"); print(io, string(UInt8(c); base=16, pad=2))
        else              print(io, c)
        end
    end
    print(io, '"')
end

json_number(io::IO, x::Float64) = isfinite(x) ? print(io, x) : print(io, "null")
json_number(io::IO, x::Int) = print(io, x)
json_bool(io::IO, b::Bool) = print(io, b ? "true" : "false")

function json_array(io::IO, v::AbstractVector{Float64})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        json_number(io, v[i])
    end
    print(io, ']')
end

function json_array(io::IO, v::AbstractVector{Int})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        print(io, v[i])
    end
    print(io, ']')
end

function json_array(io::IO, v::AbstractVector{String})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        json_string(io, v[i])
    end
    print(io, ']')
end

# `"key":` with a leading comma unless it is the first member.
function json_key(io::IO, k::String, first::Bool=false)
    first || print(io, ',')
    json_string(io, k)
    print(io, ':')
end

function json_head(io::IO, o::Options, source::String)
    print(io, '{')
    json_key(io, "schema_version", true); print(io, SCHEMA_VERSION)
    json_key(io, "sigmatau"); json_string(io, SIGMATAU_VERSION)
    json_key(io, "command"); json_string(io, o.command)
    json_key(io, "source"); json_string(io, source)
end

# The source is the one a result of the command names, as far as it is known without the record.
function write_error_json(io::IO, o::Options, code::Int, msg::String)
    json_head(io, o, o.command == "model" ? "model" : !isempty(o.input) ? o.input : o.synth > 0 ? "synthetic" : "")
    json_key(io, "error"); json_string(io, msg)
    json_key(io, "exit_code"); print(io, code)
    print(io, "}\n")
end

# Samples of x to plot with at most n points: x itself when short enough,
# otherwise its first and last sample and, in each of (n - 2) ÷ 3 equal
# stretches, the lowest and highest sample and the first gap, in time order, so
# the plot spans the record and no outlier or gap falls between drawn points.
function series_index(x::AbstractVector, n::Int)
    N = length(x)
    N <= n && return collect(1:N)
    nb = max(1, (n - 2) ÷ 3)
    idx = [1]
    for b in 1:nb
        imin = imax = inan = 0
        for i in ((b - 1) * N ÷ nb + 1):(b * N ÷ nb)
            v = x[i]
            if isnan(v)
                inan == 0 && (inan = i)
            else
                (imin == 0 || v < x[imin]) && (imin = i)
                (imax == 0 || v > x[imax]) && (imax = i)
            end
        end
        for i in sort!([imin, imax, inan])
            i > 0 && idx[end] != i && push!(idx, i)
        end
    end
    idx[end] == N || push!(idx, N)
    return idx
end

# The record and its summary; `series` adds the phase thinned for plotting
# (`series_index`).
function json_record(io::IO, o::Options, data::PhaseData)
    s = summarize(data)
    tau0 = Float64(data.tau0)
    print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "n"); print(io, s.n)
    json_key(io, "tau0"); json_number(io, tau0)
    json_key(io, "span_s"); json_number(io, Float64(s.span))
    json_key(io, "start"); json_string(io, isnan(data.start) ? "" : SigmaTau._iso_text(data.start))
    json_key(io, "gaps"); print(io, s.gaps)
    json_key(io, "outliers"); print(io, length(s.outliers))
    keep = s.outliers[1:min(length(s.outliers), 500)]
    json_key(io, "outlier_t"); json_array(io, Float64[(i - 1) * tau0 for i in keep])
    json_key(io, "outlier_x"); json_array(io, Float64[Float64(data.x[i]) for i in keep])
    json_key(io, "drift_per_day"); json_number(io, Float64(s.drift_per_day))
    json_key(io, "removed_offset"); json_number(io, o.removed_offset)
    json_key(io, "removed_drift_per_day"); json_number(io, o.removed_drift)
    if o.series > 0
        idx = series_index(data.x, o.series)
        json_key(io, "series")
        print(io, '{')
        json_key(io, "t", true); json_array(io, Float64[(i - 1) * tau0 for i in idx])
        json_key(io, "x"); json_array(io, Float64[Float64(data.x[i]) for i in idx])
        print(io, '}')
    end
    print(io, '}')
end

function json_preprocessing(io::IO, o::Options)
    json_key(io, "detrend"); json_string(io, o.detrend)
    json_key(io, "remove_drift"); json_string(io, o.remove_drift)
    json_key(io, "fillgaps"); json_bool(io, o.fillgaps)
    json_key(io, "remove_outliers"); json_bool(io, o.remove_outliers)
    json_key(io, "reference_col"); print(io, o.reference_col)
    json_key(io, "range"); print(io, '[')
    json_number(io, o.range_start); print(io, ','); json_number(io, o.range_end); print(io, ']')
end

function json_result(io::IO, r::StabilityResult, label::String, axis::String)
    has_ci = !isempty(r.ci)
    has_edf = !isempty(r.edf)
    has_nt = !isempty(r.noise_type)
    print(io, '{')
    json_key(io, "deviation", true); json_string(io, String(r.deviation_type))
    json_key(io, "label"); json_string(io, label)
    json_key(io, "axis"); json_string(io, axis)
    json_key(io, "points"); print(io, '[')
    for i in eachindex(r.tau)
        i > 1 && print(io, ',')
        print(io, '{')
        json_key(io, "tau", true); json_number(io, r.tau[i])
        json_key(io, "dev"); json_number(io, r.dev[i])
        json_key(io, "ci_lower"); has_ci ? json_number(io, r.ci[i].lo) : print(io, "null")
        json_key(io, "ci_upper"); has_ci ? json_number(io, r.ci[i].hi) : print(io, "null")
        json_key(io, "edf"); has_edf ? json_number(io, r.edf[i]) : print(io, "null")
        json_key(io, "neff"); print(io, r.neff[i])
        json_key(io, "noise_type"); has_nt ? json_string(io, String(r.noise_type[i])) : print(io, "null")
        print(io, '}')
    end
    print(io, "]}")
end

# The power-law model on MHDEV: fitted to every α from 2 to −2, or with
# --component built by hand one noise type at a time. Written: h_α, the model
# curve of each result it covers, the model ADEV at round averaging times
# within half the record, and the clock model it gives. With components, also
# each component's τ span and the MHDEV left after them. A record or a
# component the fit cannot take answers {"error": …}.
const SPEC_TAUS = [1.0, 10.0, 100.0, 1000.0, 1.0e4, 86400.0]

function json_powerlaw(io::IO, o::Options, data::PhaseData, results::Vector{StabilityResult})
    hand = !isempty(o.components)
    nf = NoiseFit[]
    h = try
        if hand
            f = noise_fit(mhdev(data, SigmaTau._default_m_values(length(data.x), :mhdev)))
            for (α, a, b) in o.components
                k = findall(t -> !(t < a) && !(t > b), f.result.tau)   # NaN bounds are open
                isempty(k) && throw(ArgumentError("--component: no MHDEV averaging time in that span"))
                add_component!(f, α, first(k):last(k))
            end
            push!(nf, f)
            Dict{Int,Float64}(c.alpha => c.h for c in f.components)
        else
            pf = powerlaw_fit(data)
            Dict{Int,Float64}(zip(pf.alphas, pf.h))
        end
    catch err
        err isa ArgumentError || rethrow()
        print(io, '{'); json_key(io, "error", true); json_string(io, error_message(err)); print(io, '}')
        return nothing
    end
    τ₀ = data.tau0
    αs = reverse!(sort!(collect(keys(h))))   # sort!(…; rev) is not trim-safe
    print(io, '{')
    json_key(io, "deviation", true); json_string(io, "mhdev")
    json_key(io, "alphas"); json_array(io, αs)
    json_key(io, "h"); json_array(io, Float64[h[α] for α in αs])
    json_key(io, "noise_types"); json_array(io, String[String(SigmaTau._nnls_symbol(a)) for a in αs])
    json_key(io, "mode"); json_string(io, hand ? "components" : "fit")
    if hand
        f = nf[1]
        json_key(io, "components"); print(io, '[')
        for (i, c) in enumerate(f.components)
            i > 1 && print(io, ',')
            print(io, '{')
            json_key(io, "alpha", true); print(io, c.alpha)
            json_key(io, "h"); json_number(io, c.h)
            json_key(io, "tau_from"); json_number(io, f.result.tau[first(c.range)])
            json_key(io, "tau_to"); json_number(io, f.result.tau[last(c.range)])
            print(io, '}')
        end
        print(io, ']')
        json_key(io, "residual"); print(io, '{')
        json_key(io, "tau", true); json_array(io, f.result.tau)
        json_key(io, "dev"); json_array(io, Float64[v > 0 ? sqrt(v) : 0.0 for v in noise_residual(f)])
        print(io, '}')
    end
    json_key(io, "curves"); print(io, '[')
    k = 0
    for r in results
        d = r.deviation_type
        base = d === :tdev ? :mdev : d
        base in SigmaTau.POWERLAW_DEVIATIONS || continue
        model = powerlaw_sigma(h, r.tau; deviation=base, tau0=τ₀)
        d === :tdev && (model .*= r.tau ./ sqrt(3.0))
        (k += 1) > 1 && print(io, ',')
        print(io, '{')
        json_key(io, "deviation", true); json_string(io, String(d))
        json_key(io, "tau"); json_array(io, r.tau)
        json_key(io, "dev"); json_array(io, model)
        print(io, '}')
    end
    print(io, ']')
    half = (length(data.x) - 1) * τ₀ / 2
    taus = unique(Float64[max(1, round(Int, t / τ₀)) * τ₀ for t in SPEC_TAUS])
    filter!(t -> t <= half, taus)
    json_key(io, "sigma_at"); print(io, '{')
    json_key(io, "tau", true); json_array(io, taus)
    json_key(io, "adev"); json_array(io, powerlaw_sigma(h, taus; deviation=:adev, tau0=τ₀))
    print(io, '}')
    m = clock_model(h; tau0=τ₀)
    json_key(io, "clock_model"); print(io, '{')
    json_key(io, "R", true); json_number(io, m.R)
    json_key(io, "q1"); json_number(io, m.q1)
    json_key(io, "q2"); json_number(io, m.q2)
    json_key(io, "q3"); json_number(io, m.q3)
    print(io, "}}")
    return nothing
end

# ITU-T wander masks (one, or every one for "all"): each one's MTIE and TDEV
# limits on a log τ grid (null where the recommendation does not reach), and
# its check of each MTIE / TDEV result.
const MASK_TAUS = [10.0^(k / 20) for k in -20:120]

function json_masks(io::IO, name::String, results::Vector{StabilityResult})
    names = name == "all" ? collect(Symbol, SigmaTau.TELECOM_MASKS) : [Symbol(name)]
    print(io, '[')
    for (j, m) in enumerate(names)
        j > 1 && print(io, ',')
        print(io, '{')
        json_key(io, "name", true); json_string(io, String(m))
        json_key(io, "title"); json_string(io, SigmaTau._MASKS[m][1])
        json_key(io, "tau"); json_array(io, MASK_TAUS)
        json_key(io, "mtie"); json_array(io, telecom_mask(m, :mtie, MASK_TAUS))
        json_key(io, "tdev"); json_array(io, telecom_mask(m, :tdev, MASK_TAUS))
        json_key(io, "checks"); print(io, '[')
        k = 0
        for r in results
            r.deviation_type in (:mtie, :tdev) || continue
            c = mask_check(r, m)
            (k += 1) > 1 && print(io, ',')
            print(io, '{')
            json_key(io, "deviation", true); json_string(io, String(r.deviation_type))
            json_key(io, "limit"); json_array(io, c.limit)
            json_key(io, "pass"); json_bool(io, c.pass)
            json_key(io, "margin"); json_number(io, c.margin)
            print(io, '}')
        end
        print(io, "]}")
    end
    print(io, ']')
    return nothing
end

function json_spectrum(io::IO, sy::SpectralResult, sx::SpectralResult)
    n = min(length(sy.freq), length(sx.freq))
    print(io, '{')
    json_key(io, "nperseg", true); print(io, sy.nperseg)
    json_key(io, "noverlap"); print(io, sy.noverlap)
    json_key(io, "f"); json_array(io, view(sx.freq, 2:n))
    json_key(io, "sy"); json_array(io, view(sy.psd, 2:n))
    json_key(io, "sx"); json_array(io, view(sx.psd, 2:n))
    print(io, '}')
end

function json_dynamic(io::IO, r::DynamicStabilityResult, step::Int)
    print(io, '{')
    json_key(io, "deviation", true); json_string(io, String(r.deviation_type))
    json_key(io, "window"); print(io, r.window)
    json_key(io, "step"); print(io, step)
    json_key(io, "t"); json_array(io, r.t)
    json_key(io, "tau"); json_array(io, r.tau)
    json_key(io, "dev"); print(io, '[')
    for i in 1:size(r.dev, 1)
        i > 1 && print(io, ',')
        json_array(io, view(r.dev, i, :))
    end
    print(io, "]}")
end

# ── CSV ───────────────────────────────────────────────────────────────────

function csv_field(io::IO, s::String)
    if any(c -> c == ',' || c == '"' || c == '\n' || c == '\r', s)
        print(io, '"'); print(io, replace(s, "\"" => "\"\"")); print(io, '"')
    else
        print(io, s)
    end
end

# Blank for NaN or ±Inf, never the text "NaN".
csv_number(io::IO, x::Float64) = isfinite(x) && print(io, x)

const CSV_ANALYZE = "source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type"

function csv_result(io::IO, source::String, q::String, r::StabilityResult)
    has_ci = !isempty(r.ci)
    has_edf = !isempty(r.edf)
    has_nt = !isempty(r.noise_type)
    for i in eachindex(r.tau)
        csv_field(io, source); print(io, ','); print(io, q); print(io, ',')
        print(io, String(r.deviation_type)); print(io, ',')
        csv_number(io, r.tau[i]); print(io, ',')
        csv_number(io, r.dev[i]); print(io, ',')
        has_ci && csv_number(io, r.ci[i].lo); print(io, ',')
        has_ci && csv_number(io, r.ci[i].hi); print(io, ',')
        has_edf && csv_number(io, r.edf[i]); print(io, ',')
        print(io, r.neff[i]); print(io, ',')
        has_nt && print(io, String(r.noise_type[i]))
        print(io, '\n')
    end
end
