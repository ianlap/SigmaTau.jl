# The analysis commands. Each computes first, then writes CSV or JSON.

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

function run_command(io::IO, o::Options, data::PhaseData, t_start::UInt64, args::Vector{String})
    c = o.command
    c == "report"   && return cmd_report(io, o, data, args)
    c == "analyze"  && return cmd_analyze(io, o, data, t_start)
    c == "check"    && return cmd_check(io, o, data, t_start)
    c == "spectrum" && return cmd_spectrum(io, o, data, t_start)
    c == "dynamic"  && return cmd_dynamic(io, o, data, t_start)
    c == "cross"    && return cmd_cross(io, o, data, t_start)
    c == "drift"    && return cmd_drift(io, o, data, t_start)
    c == "holdover" && return cmd_holdover(io, o, data, t_start)
    c == "steer"    && return cmd_steer(io, o, data, t_start)
    c == "timescale" && return cmd_timescale(io, o, data, t_start)
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

# One method instance per type of `bias` and `noise`: under --trim=safe the
# keyword call needs concrete ones.
run_deviation(d::Symbol, data, ms::Vector{Int}, o::Options, bias::Union{Bool,Symbol}, noise) =
    deviation(d, data, ms; ci=o.ci, confidence=o.confidence, correct_bias=bias, noise=noise)

function analysis_results(o::Options, data::PhaseData)
    need_samples(data, 8, "the analysis")
    o.tau_min > o.tau_max && throw(usage_error("--tau-min is above --tau-max"))
    N = length(data.x)
    mode = tau_mode(o.grid)
    results = StabilityResult[]
    fit = o.noise == "fit" ? powerlaw_fit(data) : nothing
    for d in o.devs
        ms = isempty(o.m) ? tau_values(mode, N, d) : o.m
        ms = Int[m for m in ms if !(m * data.tau0 < o.tau_min * (1 - 1e-9)) && !(m * data.tau0 > o.tau_max * (1 + 1e-9))]
        isempty(ms) && continue
        push!(results, fit === nothing ? run_deviation(d, data, ms, o, correct_bias(o), Symbol(o.noise)) :
                                         run_deviation(d, data, ms, o, correct_bias(o), fit))
    end
    isempty(results) && !(isnan(o.tau_min) && isnan(o.tau_max)) &&
        throw(CLIError(EXIT_ANALYSIS, "no averaging time of the grid lies within --tau-min and --tau-max"))
    return results
end

function cmd_analyze(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    results = analysis_results(o, data)
    N = length(data.x)
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
    json_key(io, "tau_min"); json_number(io, o.tau_min)
    json_key(io, "tau_max"); json_number(io, o.tau_max)
    json_key(io, "ci"); json_bool(io, o.ci)
    json_key(io, "confidence"); o.ci ? json_number(io, o.confidence) : print(io, "null")
    json_key(io, "bias_correction"); json_string(io, o.bias)
    json_key(io, "noise"); json_string(io, o.noise)
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
        if any(isnan, data.x)
            # The estimators skip gaps; the spectrum cannot, so it is left out.
            print(io, '{'); json_key(io, "error", true); json_string(io, "the spectrum needs a record without gaps"); print(io, '}')
        else
            sy, sx = spectra(o, data)
            json_spectrum(io, sy, sx)
        end
    end
    o.fit && (json_key(io, "powerlaw"); json_powerlaw(io, o, data, results))
    isempty(o.mask) || (json_key(io, "masks"); json_masks(io, o.mask, results))
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
    # Window centres in the whole record's time, not the range's.
    r.t .+= range_offset(o, Float64(data.tau0))
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

# ── drift ─────────────────────────────────────────────────────────────────

const CSV_DRIFT = "source,quantity,model,a,a_err,b,b_err,c,c_err,drift_per_day,drift_per_day_err,rms,n"

function cmd_drift(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 3, "a drift fit")   # drift_fit checks each model's own minimum
    f = drift_fit(data; model=Symbol(o.model))
    src = source_name(o, data)
    n = count(isfinite, f.y)
    if o.format != "json"
        print(io, CSV_DRIFT); print(io, '\n')
        csv_field(io, src); print(io, ','); print(io, quantity(o)); print(io, ','); print(io, o.model)
        for k in 1:3
            print(io, ','); k <= length(f.coef) && csv_number(io, f.coef[k])
            print(io, ','); k <= length(f.coef) && csv_number(io, f.stderr[k])
        end
        print(io, ','); csv_number(io, f.rate); print(io, ','); csv_number(io, f.rate_stderr)
        print(io, ','); csv_number(io, f.rms); print(io, ','); print(io, n); print(io, '\n')
        return nothing
    end
    json_head(io, o, src)
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_key(io, "model"); json_string(io, o.model)
    json_preprocessing(io, o); print(io, '}')
    json_key(io, "drift"); print(io, '{')
    json_key(io, "model", true); json_string(io, o.model)
    json_key(io, "coef"); json_array(io, f.coef)
    json_key(io, "stderr"); json_array(io, f.stderr)
    json_key(io, "drift_per_day"); json_number(io, f.rate)
    json_key(io, "drift_per_day_err"); json_number(io, f.rate_stderr)
    json_key(io, "rms"); json_number(io, f.rms)
    json_key(io, "n"); print(io, n)
    if o.series > 0
        idx = collect(1:max(1, cld(length(f.t), o.series)):length(f.t))
        json_key(io, "series"); print(io, '{')
        json_key(io, "t", true); json_array(io, f.t[idx])
        json_key(io, "y"); json_array(io, f.y[idx])
        json_key(io, "fit"); json_array(io, f.fit[idx])
        print(io, '}')
    end
    print(io, '}')
    return elapsed(io, t_start)
end

# ── holdover ──────────────────────────────────────────────────────────────

const CSV_HOLDOVER = "source,train,tau,rms,quantile,worst,mean,n,limit,time_to_limit,claimed"

# The error the power-law model fitted to the record predicts for the same
# test, and each noise type's share of its variance.
function json_holdover_fit(io::IO, data::PhaseData, h::HoldoverResult)
    p = try
        pf = powerlaw_fit(data)
        hα, τ₀, limit = Dict{Int,Float64}(zip(pf.alphas, pf.h)), Float64(data.tau0), isnan(h.limit) ? nothing : h.limit
        h.grow > 0 ? holdover(hα; tau0=τ₀, grow=h.grow, tau=h.tau, drift=h.drift, anchor=h.anchor, limit, confidence=h.confidence) :
                     holdover(hα; tau0=τ₀, train=h.train, tau=h.tau, drift=h.drift, anchor=h.anchor, limit, confidence=h.confidence)
    catch err
        err isa ArgumentError || rethrow()
        print(io, '{'); json_key(io, "error", true); json_string(io, error_message(err)); print(io, '}')
        return nothing
    end
    print(io, '{')
    json_key(io, "alphas", true); json_array(io, p.alphas)
    json_key(io, "noise_types"); json_array(io, String[String(SigmaTau._nnls_symbol(a)) for a in p.alphas])
    json_key(io, "sigma"); json_array(io, p.sigma)
    json_key(io, "quantile"); json_array(io, p.quantile)
    json_key(io, "budget"); print(io, '[')
    for j in eachindex(p.alphas)
        j > 1 && print(io, ',')
        json_array(io, p.budget[:, j])
    end
    print(io, ']')
    json_key(io, "time_to_limit"); json_number(io, p.time_to_limit)
    print(io, '}')
    return nothing
end

# The error of the same extrapolation made from τ-averages of the phase, two
# for a line and three for a parabola: √6 TDEV and √20 HTDEV by definition.
# NaN past the τ the record gives them.
function holdover_reference(data::PhaseData, h::HoldoverResult)
    ms = Int[round(Int, t / data.tau0) for t in h.tau]
    return h.drift ? sqrt(20.0) .* htdev(data, ms; ci=false).dev : sqrt(6.0) .* tdev(data, ms; ci=false).dev
end

# The clock model of a Kalman predictor: the power-law fit to the first
# training stretch alone, so no trial is predicted with noise levels learned
# from the readings it is compared with. Flicker the fit finds is kept, as
# Markov states out to the record's span: without it the filter of a flicker
# FM clock claims a tenth of the error it makes.
function holdover_model(data::PhaseData, train::Float64, drift::Bool)
    τ₀ = Float64(data.tau0)
    b = min(SigmaTau._whole_samples(train, τ₀, "train"), length(data.x))
    pf = powerlaw_fit(PhaseData(data.x[1:b], data.tau0))
    return clock_model(Dict{Int,Float64}(zip(pf.alphas, pf.h)); tau0=τ₀, drift, flicker=true, tau_max=length(data.x) * τ₀)
end

function cmd_holdover(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 16, "a holdover test")
    τ₀ = Float64(data.tau0)
    train = isnan(o.train) ? (length(data.x) ÷ 4) * τ₀ : o.train
    limit = isnan(o.limit) ? nothing : o.limit
    tau = isempty(o.taus) ? nothing : o.taus
    model = o.kalman ? holdover_model(data, train, o.fit_drift) : nothing
    # One call per predictor and kind of stretch: the trimmed build resolves each on its own.
    h = model !== nothing ? holdover(data, model; train, tau, limit, confidence=o.confidence) :
        o.grow > 0 ? holdover(data; grow=o.grow, tau, drift=o.fit_drift, anchor=o.anchor, limit, confidence=o.confidence) :
                     holdover(data; train, tau, drift=o.fit_drift, anchor=o.anchor, limit, confidence=o.confidence)
    src = source_name(o, data)
    if o.format != "json"
        print(io, CSV_HOLDOVER); print(io, '\n')
        for i in eachindex(h.tau)
            csv_field(io, src); print(io, ','); csv_number(io, h.grow > 0 ? h.grow * h.tau[i] : h.train); print(io, ',')
            csv_number(io, h.tau[i]); print(io, ','); csv_number(io, h.rms[i]); print(io, ',')
            csv_number(io, h.quantile[i]); print(io, ','); csv_number(io, h.worst[i]); print(io, ',')
            csv_number(io, h.mean[i]); print(io, ','); print(io, h.n[i]); print(io, ',')
            csv_number(io, h.limit); print(io, ','); csv_number(io, h.time_to_limit); print(io, ',')
            model === nothing || csv_number(io, h.claimed[i]); print(io, '\n')
        end
        return nothing
    end
    json_head(io, o, src)
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_preprocessing(io, o); print(io, '}')
    json_key(io, "holdover"); print(io, '{')
    json_key(io, "train", true); json_number(io, h.train)
    json_key(io, "grow"); print(io, h.grow)
    json_key(io, "drift"); json_bool(io, h.drift)
    json_key(io, "anchor"); json_bool(io, h.anchor)
    json_key(io, "confidence"); json_number(io, h.confidence)
    json_key(io, "tau"); json_array(io, h.tau)
    json_key(io, "rms"); json_array(io, h.rms)
    json_key(io, "quantile"); json_array(io, h.quantile)
    json_key(io, "worst"); json_array(io, h.worst)
    json_key(io, "mean"); json_array(io, h.mean)
    json_key(io, "n"); json_array(io, h.n)
    json_key(io, "reference"); json_array(io, holdover_reference(data, h))
    json_key(io, "limit"); json_number(io, h.limit)
    json_key(io, "time_to_limit"); json_number(io, h.time_to_limit)
    o.fit && (json_key(io, "predicted"); json_holdover_fit(io, data, h))
    if model !== nothing
        json_key(io, "claimed"); json_array(io, h.claimed)
        json_key(io, "model"); json_clock_model(io, model)
    end
    print(io, '}')
    return elapsed(io, t_start)
end

# ── model ─────────────────────────────────────────────────────────────────

# σ(τ) the --noise-model clocks predict (MDEV, ADEV, MHDEV) at the τ of an N = --synth record.
function run_model(io::IO, o::Options)
    isempty(o.models) && throw(usage_error("model: give a clock with --noise-model R:Q1:Q2:Q3"))
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    ms = tau_values(tau_mode(o.grid), o.synth, :mdev)
    tau = tau0 .* ms
    json_head(io, o, "model")
    json_key(io, "model"); print(io, '{')
    json_key(io, "tau0", true); json_number(io, tau0)
    json_key(io, "n"); print(io, o.synth)
    json_key(io, "tau"); json_array(io, tau)
    json_key(io, "clocks"); print(io, '[')
    for (i, (R, q1, q2, q3)) in enumerate(o.models)
        i > 1 && print(io, ',')
        m = ClockModel(; tau0, R, q1, q2, q3)
        print(io, '{')
        # Random-run FM (q3) diverges in ADEV and MDEV; only the Hadamard family converges.
        for (j, d) in enumerate((:mdev, :adev, :mhdev))
            json_key(io, String(d), j == 1)
            q3 > 0 && d !== :mhdev ? print(io, "null") : json_array(io, powerlaw_sigma(m, tau; deviation=d))
        end
        print(io, '}')
    end
    print(io, "]}}\n")
    return nothing
end

# ── steer ─────────────────────────────────────────────────────────────────

const CSV_STEER = "source,tau,free,steered,predicted"

# The steady state the clock model predicts for the loop, or why there is none.
function steer_prediction(model::ClockModel, g::SteeringGains, tau::Vector{Float64})
    return try
        steer(model, g; tau, readings=true)
    catch err
        err isa ArgumentError || rethrow()
        error_message(err)
    end
end

function rms_finite(v::AbstractVector{Float64})
    s, n = 0.0, 0
    for x in v
        isnan(x) && continue
        s += abs2(x); n += 1
    end
    return n == 0 ? NaN : sqrt(s / n)
end

function json_clock_model(io::IO, m::ClockModel)
    print(io, '{')
    json_key(io, "R", true); json_number(io, m.R)
    json_key(io, "q1"); json_number(io, m.q1)
    json_key(io, "q2"); json_number(io, m.q2)
    json_key(io, "q3"); json_number(io, m.q3)
    json_key(io, "flicker_fm"); json_number(io, m.flicker_fm)
    json_key(io, "flicker_pm"); json_number(io, m.flicker_pm)
    print(io, '}')
end

# The record is the free-running clock against its reference. The clock model
# is fitted to the first training stretch alone, as for the Kalman holdover
# test; the loop runs over the whole record and every result is taken from
# the end of that stretch on.
function cmd_steer(io::IO, o::Options, data::PhaseData, t_start::UInt64)
    need_samples(data, 64, "a steering replay")
    τ₀, N = Float64(data.tau0), length(data.x)
    train = isnan(o.train) ? (N ÷ 4) * τ₀ : o.train
    b = SigmaTau._whole_samples(train, τ₀, "train")
    16 <= b <= N - 16 || throw(CLIError(EXIT_ANALYSIS, "steer: the training stretch must leave 16 samples on either side of its end"))
    interval = isnan(o.interval) ? τ₀ : o.interval
    g = isempty(o.cost) ? steering_gains(o.time_constant, interval; integral=o.integral) :
        lqg_gains(interval; phase=o.cost[1], frequency=o.cost[2], steer=o.cost[3])
    model = holdover_model(data, train, o.fit_drift)
    r = steer(data, model, g)
    free = PhaseData(Float64[Float64(x) for x in view(data.x, b+1:N)], τ₀)
    steered = PhaseData(r.steered.x[b+1:N], τ₀)
    ms = isempty(o.m) ? tau_values(tau_mode(o.grid), N - b, :adev) : o.m
    σfree, σsteered = adev(free, ms; ci=false), adev(steered, ms; ci=false)
    p = steer_prediction(model, g, Float64[m * τ₀ for m in ms])
    src = source_name(o, data)
    if o.format != "json"
        print(io, CSV_STEER); print(io, '\n')
        for i in eachindex(ms)
            csv_field(io, src); print(io, ','); csv_number(io, σfree.tau[i]); print(io, ',')
            csv_number(io, σfree.dev[i]); print(io, ','); csv_number(io, σsteered.dev[i]); print(io, ',')
            p isa SteeringPrediction && csv_number(io, p.adev[i]); print(io, '\n')
        end
        return nothing
    end
    poles = steering_poles(g)
    n = round(Int, interval / τ₀)
    json_head(io, o, src)
    json_key(io, "record"); json_record(io, o, data)
    json_key(io, "analysis"); print(io, '{')
    json_key(io, "quantity", true); json_string(io, quantity(o))
    json_preprocessing(io, o); print(io, '}')
    json_key(io, "steer"); print(io, '{')
    json_key(io, "train", true); json_number(io, train)
    json_key(io, "interval"); json_number(io, interval)
    json_key(io, "gains"); print(io, '{')
    json_key(io, "phase", true); json_number(io, g.phase)
    json_key(io, "frequency"); json_number(io, g.frequency)
    json_key(io, "integral"); json_number(io, g.integral)
    print(io, '}')
    json_key(io, "time_constants"); json_array(io, poles.time_constants)
    json_key(io, "frequencies"); json_array(io, poles.frequencies)
    json_key(io, "drift"); json_bool(io, model.drift)
    json_key(io, "model"); json_clock_model(io, model)
    json_key(io, "tau"); json_array(io, σfree.tau)
    json_key(io, "free"); json_array(io, σfree.dev)
    json_key(io, "steered"); json_array(io, σsteered.dev)
    json_key(io, "phase_rms"); json_number(io, rms_finite(steered.x))
    json_key(io, "steers_rms"); json_number(io, rms_finite(r.steers[(b ÷ n + 1) * n:n:N]))
    json_key(io, "predicted"); print(io, '{')
    if p isa SteeringPrediction
        json_key(io, "phase", true); json_number(io, p.phase)
        json_key(io, "frequency"); json_number(io, p.frequency)
        json_key(io, "steers"); json_number(io, p.steers)
        json_key(io, "adev"); json_array(io, p.adev)
    else
        json_key(io, "error", true); json_string(io, p)
    end
    print(io, '}')
    if o.series > 0
        idx = series_index(r.steered.x, o.series)
        json_key(io, "series"); print(io, '{')
        json_key(io, "t", true); json_array(io, Float64[(i - 1) * τ₀ for i in idx])
        json_key(io, "steered"); json_array(io, r.steered.x[idx])
        json_key(io, "correction"); json_array(io, r.correction[idx])
        print(io, '}')
    end
    print(io, '}')
    return elapsed(io, t_start)
end

# ── timescale ─────────────────────────────────────────────────────────────

# White + random-walk FM clocks (four unless --clocks says otherwise), each a
# factor of two above the last in level, read against a perfect reference.
function synth_ensemble(o::Options)
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    ms = if isempty(o.models)
        [ClockModel(; tau0, q1=1e-24 * 4.0^j, q2=1e-30 * 4.0^j) for j in 0:(o.clocks > 0 ? o.clocks : 4) - 1]
    else
        [ClockModel(; tau0, R, q1, q2, q3) for (R, q1, q2, q3) in o.models]
    end
    cc, _ = simulate_ensemble(ms, o.synth, reference=0, rng=Random.Xoshiro(o.seed))
    return cc
end

# Each --value-cols column read like FILE, with the same options.
function read_ensemble(o::Options)
    isempty(o.value_cols) && throw(usage_error("timescale needs --value-cols with two or more columns, or --synth"))
    col, offset, drift = o.value_col, o.removed_offset, o.removed_drift
    records = PhaseData{Float64}[]
    try
        for c in o.value_cols
            o.value_col = c
            push!(records, load_record(o))
        end
    catch e
        e isa CLIError && rethrow()
        throw(CLIError(EXIT_INPUT, error_message(e)))   # read after the input phase
    finally
        o.value_col, o.removed_offset, o.removed_drift = col, offset, drift
    end
    allequal(length(r.x) for r in records) && allequal(r.tau0 for r in records) ||
        throw(CLIError(EXIT_INPUT, "the clock columns differ in length or τ₀"))
    return ClockComparisons(records; names=["column " * string(c) for c in o.value_cols])
end

# The Kalman scales fit a model per clock (with --flicker, flicker states
# included) and read the scale at --pivot. Returns the scale and the models.
function run_timescale(cc::ClockComparisons, o::Options, kalman::Bool)
    alg = Symbol(o.algorithm)
    o.flicker && !(alg === :greenhall || alg === :jones_tryon) && throw(usage_error("--flicker needs --algorithm greenhall or jones_tryon"))
    if !kalman
        isempty(o.steps) || alg === :at2 || throw(usage_error("--steps needs --algorithm at2"))
        isempty(o.steps) && return timescale(cc; algorithm=alg, wmax=o.wmax), ClockModel[]
        return timescale(cc; algorithm=alg, wmax=o.wmax, steps=Symbol(o.steps)), ClockModel[]
    end
    k = size(cc.x, 2)
    k >= 3 || throw(usage_error("--algorithm $(o.algorithm) fits a model per clock and needs three or more clocks"))
    1 <= o.pivot <= k || throw(usage_error("--pivot must be between 1 and $k"))
    models = clock_models(cc; flicker=o.flicker)
    return timescale(cc; algorithm=alg, models, pivot=o.pivot, wmax=o.wmax), models
end

function cmd_timescale(io::IO, o::Options, ::PhaseData, t_start::UInt64)
    cc = o.synth > 0 ? synth_ensemble(o) : read_ensemble(o)
    N, k = size(cc.x)
    N >= 16 || throw(CLIError(EXIT_ANALYSIS, "a timescale needs at least 16 samples"))
    kalman = o.algorithm in ("greenhall", "jones_tryon", "stein")
    r, models = run_timescale(cc, o, kalman)
    late = max(1, N ÷ 2):N
    wmean = [sum(r.weights[late, i]) / length(late) for i in 1:k]
    rms = Float64[]
    for i in 1:k
        e = filter(isfinite, r.error[late, i])
        push!(rms, isempty(e) ? NaN : sqrt(sum(abs2, e) / length(e)))
    end
    if o.format != "json"
        print(io, "clock,weight,weight_mean,error_rms\n")
        for i in 1:k
            csv_field(io, cc.names[i]); print(io, ',')
            csv_number(io, r.weights[end, i]); print(io, ',')
            csv_number(io, wmean[i]); print(io, ',')
            csv_number(io, rms[i]); print(io, '\n')
        end
        return nothing
    end
    ms = tau_values(tau_mode(o.grid), N, :adev)
    json_head(io, o, o.synth > 0 ? "synthetic" : o.input)
    json_key(io, "timescale"); print(io, '{')
    json_key(io, "algorithm", true); json_string(io, o.algorithm)
    json_key(io, "tau0"); json_number(io, cc.tau0)
    json_key(io, "n"); print(io, N)
    json_key(io, "reference"); json_string(io, cc.reference)
    json_key(io, "wmax"); json_number(io, o.wmax)
    kalman && (json_key(io, "pivot"); json_number(io, o.pivot))
    isempty(o.steps) || (json_key(io, "steps"); json_string(io, o.steps))
    if o.flicker   # Markov states each clock's model carries
        json_key(io, "flicker"); print(io, '{')
        json_key(io, "fm", true); json_array(io, Int[length(SigmaTau._fm_taus(m)) for m in models])
        json_key(io, "pm"); json_array(io, Int[length(SigmaTau._pm_taus(m)) for m in models])
        print(io, '}')
    end
    json_key(io, "clocks"); json_array(io, cc.names)
    json_key(io, "weight"); json_array(io, r.weights[end, :])
    json_key(io, "weight_mean"); json_array(io, wmean)
    json_key(io, "error_rms"); json_array(io, rms)
    # ADEV of the scale against the reference, and of each clock against the scale.
    json_key(io, "adev"); print(io, '{')
    s = adev(timescale_phase(r), ms; ci=false)
    json_key(io, "tau", true); json_array(io, s.tau)
    json_key(io, "scale"); json_array(io, s.dev)
    json_key(io, "clocks"); print(io, '[')
    for i in 1:k
        i > 1 && print(io, ',')
        json_array(io, adev(PhaseData(r.x[:, i], cc.tau0), ms; ci=false).dev)
    end
    print(io, "]}")
    if o.series > 0
        idx = collect(1:max(1, cld(N, o.series)):N)
        json_key(io, "series"); print(io, '{')
        json_key(io, "t", true); json_array(io, r.t[idx])
        json_key(io, "scale"); json_array(io, timescale_phase(r).x[idx])
        json_key(io, "weights"); print(io, '[')
        for i in 1:k
            i > 1 && print(io, ',')
            json_array(io, r.weights[idx, i])
        end
        print(io, "]}")
    end
    print(io, '}')
    return elapsed(io, t_start)
end
