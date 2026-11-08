# The four analysis commands. Each computes first, then writes CSV or JSON.

source_name(o::Options, data::PhaseData) = o.synth > 0 ? "synthetic" : o.input

function registry_entry(name::Symbol)
    for e in estimators()
        e.name === name && return e
    end
    throw(usage_error("unknown deviation " * String(name)))
end

function need_samples(data::PhaseData, n::Int, what::String)
    N = length(data.x)
    N >= n || throw(CLIError(EXIT_ANALYSIS, what * " needs at least " * string(n) *
                                            " samples; the record has " * string(N)))
end

function elapsed(io::IO, t_start::UInt64)
    json_key(io, "elapsed_ms"); json_number(io, (time_ns() - t_start) / 1e6)
    print(io, "}\n")
end

function run_command(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    c = o.command
    c == "analyze"  && return cmd_analyze(io, o, data, t_start)
    c == "check"    && return cmd_check(io, o, data, t_start)
    c == "spectrum" && return cmd_spectrum(io, o, data, t_start)
    c == "dynamic"  && return cmd_dynamic(io, o, data, t_start)
    c == "cross"    && return cmd_cross(io, o, data, t_start)
    throw(usage_error("unknown command " * c))
end

# ── check ─────────────────────────────────────────────────────────────────

function cmd_check(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 3, "a record")
    src = source_name(o, data)
    if o.format == "json"
        json_head(io, o, src)
        json_key(io, "record"); json_record(io, o, data)
        json_key(io, "analysis"); print(io, '{')
        json_key(io, "quantity", true); json_string(io, quantity(o))
        json_preprocessing(io, o); print(io, '}')
        return elapsed(io, t_start)
    end
    s = summarize(data)
    print(io, "source,quantity,n,tau0,span_s,gaps,outliers,drift_per_day\n")
    csv_field(io, src); print(io, ','); print(io, quantity(o)); print(io, ',')
    print(io, s.n); print(io, ','); csv_number(io, Float64(s.tau0)); print(io, ',')
    csv_number(io, Float64(s.span)); print(io, ','); print(io, s.gaps); print(io, ',')
    print(io, length(s.outliers)); print(io, ','); csv_number(io, Float64(s.drift_per_day))
    print(io, '\n')
end

# ── analyze ───────────────────────────────────────────────────────────────

function cmd_analyze(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 8, "the analysis")
    N = length(data.x)
    mode = tau_mode(o.grid)
    results = StabilityResult[]
    for d in o.devs
        ms = isempty(o.m) ? tau_values(mode, N, d) : o.m
        isempty(ms) && continue
        push!(results, deviation(d, data, ms; ci=o.ci, confidence=o.confidence, correct_bias=o.bias))
    end
    src = source_name(o, data)
    if o.format != "json"
        print(io, CSV_ANALYZE); print(io, '\n')
        for r in results
            csv_result(io, src, quantity(o), r)
        end
        return nothing
    end
    json_head(io, o, src)
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "tau0"); json_number(io, Float64(data.tau0))
    json_key(io, "n"); print(io, N)
    json_key(io, "tau_mode"); json_string(io, isempty(o.m) ? o.grid : "m")
    isempty(o.m) || (json_key(io, "m"); json_array(io, o.m))
    json_key(io, "ci"); json_bool(io, o.ci)
    json_key(io, "confidence"); o.ci ? json_number(io, o.confidence) : print(io, "null")
    json_key(io, "bias_correction"); json_bool(io, o.bias)
    json_key(io, "deviations"); json_array(io, String[String(d) for d in o.devs])
    json_preprocessing(io, o)
    print(io, '}')
    json_key(io, "results"); print(io, '[')
    for (k, r) in enumerate(results)
        k > 1 && print(io, ',')
        e = registry_entry(r.deviation_type)
        json_result(io, r, e.label, String(e.axis))
    end
    print(io, ']')
    if o.spectrum
        json_key(io, "spectrum")
        sy, sx = spectra(o, data)
        json_spectrum(io, sy, sx)
    end
    return elapsed(io, t_start)
end

# ── spectrum ──────────────────────────────────────────────────────────────

function spectra(o::Options, data::PhaseData)
    need_samples(data, 16, "the spectrum")
    any(isnan, data.x) &&
        throw(CLIError(EXIT_ANALYSIS, "the spectrum needs a record without gaps; fill them with --fillgaps"))
    N = length(data.x)
    # One segment length for both spectra, from the phase record's length (the
    # library's own default for psd_y counts the N − 1 frequency samples).
    np = o.nperseg > 0 ? min(o.nperseg, N) : SigmaTau._default_nperseg(N)
    return psd_y(data; nperseg=np), psd_x(data; nperseg=np)
end

function cmd_spectrum(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    sy, sx = spectra(o, data)
    if o.format != "json"
        print(io, "f,sy,sx\n")
        for i in 2:min(length(sy.freq), length(sx.freq))
            csv_number(io, sx.freq[i]); print(io, ',')
            csv_number(io, sy.psd[i]); print(io, ',')
            csv_number(io, sx.psd[i]); print(io, '\n')
        end
        return nothing
    end
    json_head(io, o, source_name(o, data))
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "taper"); json_string(io, "hann")
    json_preprocessing(io, o); print(io, '}')
    json_key(io, "spectrum"); json_spectrum(io, sy, sx)
    return elapsed(io, t_start)
end

# ── dynamic ───────────────────────────────────────────────────────────────

function cmd_dynamic(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 16, "a dynamic deviation")
    N = length(data.x)
    w = o.window > 0 ? min(o.window, N) : SigmaTau._default_window(N)
    st = o.step > 0 ? o.step : max(w ÷ 2, 1)
    hd = o.dev == "dhdev"
    ms = tau_values(tau_mode(o.grid), w, hd ? :hdev : :adev)
    isempty(ms) && throw(CLIError(EXIT_ANALYSIS, "the window is too short for any τ"))
    r = hd ? dhdev(data, ms; window=w, step=st) : dadev(data, ms; window=w, step=st)
    if o.format != "json"
        print(io, "t,tau,dev\n")
        for i in eachindex(r.t), j in eachindex(r.tau)
            csv_number(io, r.t[i]); print(io, ',')
            csv_number(io, r.tau[j]); print(io, ',')
            csv_number(io, r.dev[i, j]); print(io, '\n')
        end
        return nothing
    end
    json_head(io, o, source_name(o, data))
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "tau_mode"); json_string(io, o.grid)
    json_preprocessing(io, o); print(io, '}')
    json_key(io, "dynamic"); json_dynamic(io, r, st)
    return elapsed(io, t_start)
end

# ── cross ─────────────────────────────────────────────────────────────────

# The second channel, read like the first.
function second_channel(o::Options)
    # The primary record's settings and drift fit stay as they were.
    first, measured, offset, drift = o.input, o.measured, o.removed_offset, o.removed_drift
    o.input, o.measured = o.with, falses(0)
    try
        return load_record(o)
    catch e
        e isa CLIError && rethrow()
        throw(CLIError(EXIT_INPUT, "--with " * o.with * ": " * error_message(e)))
    finally
        o.input, o.measured, o.removed_offset, o.removed_drift = first, measured, offset, drift
    end
end

function cmd_cross(io::IO, o::Options, a::PhaseData, t_start::UInt64)
    b = second_channel(o)
    need_samples(a, 16, "a cross analysis")
    length(a.x) == length(b.x) && a.tau0 == b.tau0 ||
        throw(CLIError(EXIT_INPUT, "the two channels differ in length or τ₀"))
    (any(isnan, a.x) || any(isnan, b.x)) &&
        throw(CLIError(EXIT_ANALYSIS, "cross analysis needs records without gaps; fill them with --fillgaps"))
    N = length(a.x)
    ms = isempty(o.m) ? tau_values(tau_mode(o.grid), N, :adev) : o.m
    r = cross_adev(a, b, ms)
    np = o.nperseg > 0 ? min(o.nperseg, N - 1) : SigmaTau._default_nperseg(N - 1)
    s = cross_psd_y(a, b; nperseg=np)
    if o.format != "json"
        print(io, "tau,covar,dev,uncertainty,dev_a,dev_b,neff\n")
        for i in eachindex(r.tau)
            csv_number(io, r.tau[i]); print(io, ',')
            csv_number(io, r.covar[i]); print(io, ',')
            csv_number(io, r.dev[i]); print(io, ',')
            csv_number(io, r.uncertainty[i]); print(io, ',')
            csv_number(io, r.dev_a[i]); print(io, ',')
            csv_number(io, r.dev_b[i]); print(io, ',')
            print(io, r.neff[i]); print(io, '\n')
        end
        return nothing
    end
    json_head(io, o, o.input)
    json_key(io, "with"); json_string(io, o.with)
    json_key(io, "record"); json_record(io, o, a)
    json_key(io, "cross_adev"); print(io, '{')
    json_key(io, "tau", true); json_array(io, r.tau)
    json_key(io, "covar"); json_array(io, r.covar)
    json_key(io, "dev"); json_array(io, r.dev)
    json_key(io, "uncertainty"); json_array(io, r.uncertainty)
    json_key(io, "dev_a"); json_array(io, r.dev_a)
    json_key(io, "dev_b"); json_array(io, r.dev_b)
    json_key(io, "neff"); json_array(io, r.neff)
    print(io, '}')
    json_key(io, "cross_spectrum"); print(io, '{')
    json_key(io, "nperseg", true); print(io, s.nperseg)
    json_key(io, "nseg"); print(io, s.nseg)
    json_key(io, "f"); json_array(io, view(s.freq, 2:length(s.freq)))
    json_key(io, "sy"); json_array(io, view(s.psd, 2:length(s.freq)))
    json_key(io, "floor"); json_array(io, view(s.floor, 2:length(s.freq)))
    json_key(io, "sy_a"); json_array(io, view(s.psd_a, 2:length(s.freq)))
    json_key(io, "sy_b"); json_array(io, view(s.psd_b, 2:length(s.freq)))
    print(io, '}')
    return elapsed(io, t_start)
end
