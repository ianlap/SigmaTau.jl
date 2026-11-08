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

function write_error_json(io::IO, o::Options, code::Int, msg::String)
    json_head(io, o, isempty(o.input) ? (o.synth > 0 ? "synthetic" : "") : o.input)
    json_key(io, "error"); json_string(io, msg)
    json_key(io, "exit_code"); print(io, code)
    print(io, "}\n")
end

# The record and its summary; `series` adds an evenly strided subset of the
# phase for plotting.
function json_record(io::IO, o::Options, data::PhaseData)
    s = summarize(data)
    tau0 = Float64(data.tau0)
    print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "n"); print(io, s.n)
    json_key(io, "tau0"); json_number(io, tau0)
    json_key(io, "span_s"); json_number(io, Float64(s.span))
    json_key(io, "gaps"); print(io, s.gaps)
    json_key(io, "outliers"); print(io, length(s.outliers))
    keep = s.outliers[1:min(length(s.outliers), 500)]
    json_key(io, "outlier_t"); json_array(io, Float64[(i - 1) * tau0 for i in keep])
    json_key(io, "outlier_x"); json_array(io, Float64[Float64(data.x[i]) for i in keep])
    json_key(io, "drift_per_day"); json_number(io, Float64(s.drift_per_day))
    if o.series > 0
        N = length(data.x)
        idx = collect(1:max(1, cld(N, o.series)):N)
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
    json_key(io, "fillgaps"); json_bool(io, o.fillgaps)
    json_key(io, "remove_outliers"); json_bool(io, o.remove_outliers)
    json_key(io, "reference_col"); print(io, o.reference_col)
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
