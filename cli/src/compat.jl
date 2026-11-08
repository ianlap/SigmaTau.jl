# Stand-ins for the library's estimator registry (estimators, deviation,
# tau_mode, summarize) and its Welch default, for SigmaTau versions that do not
# have them yet. Each is defined only when the library lacks it, so once the
# registry is in SigmaTau this file does nothing and can be deleted.

@static if !isdefined(SigmaTau, :estimators)
    const _REGISTRY = [
        (name=:adev, label="ADEV", family=:allan, axis=:sigma_y, correct_bias=false, description=""),
        (name=:mdev, label="MDEV", family=:allan, axis=:sigma_y, correct_bias=false, description=""),
        (name=:tdev, label="TDEV", family=:allan, axis=:sigma_x, correct_bias=false, description=""),
        (name=:hdev, label="HDEV", family=:allan, axis=:sigma_y, correct_bias=false, description=""),
        (name=:mhdev, label="MHDEV", family=:allan, axis=:sigma_y, correct_bias=false, description=""),
        (name=:htdev, label="HTDEV", family=:allan, axis=:sigma_x, correct_bias=false, description=""),
        (name=:totdev, label="TOTDEV", family=:total, axis=:sigma_y, correct_bias=true, description=""),
        (name=:mtotdev, label="MTOTDEV", family=:total, axis=:sigma_y, correct_bias=true, description=""),
        (name=:ttotdev, label="TTOTDEV", family=:total, axis=:sigma_x, correct_bias=true, description=""),
        (name=:htotdev, label="HTOTDEV", family=:total, axis=:sigma_y, correct_bias=true, description=""),
        (name=:mhtotdev, label="MHTOTDEV", family=:total, axis=:sigma_y, correct_bias=true, description=""),
        (name=:theo1, label="Thêo1", family=:theo, axis=:sigma_y, correct_bias=true, description=""),
        (name=:theoh, label="ThêoH", family=:theo, axis=:sigma_y, correct_bias=false, description=""),
        (name=:pdev, label="PDEV", family=:allan, axis=:sigma_y, correct_bias=false, description=""),
        (name=:mtie, label="MTIE", family=:time, axis=:sigma_x, correct_bias=false, description=""),
        (name=:tierms, label="TIE rms", family=:time, axis=:sigma_x, correct_bias=false, description=""),
    ]
    estimators() = _REGISTRY

    function deviation(name::Symbol, data::PhaseData, m::Vector{Int};
                       ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE, correct_bias::Bool=true)
        c, b = confidence, correct_bias
        name === :adev     && return adev(data, m; ci=ci, confidence=c)
        name === :mdev     && return mdev(data, m; ci=ci, confidence=c)
        name === :tdev     && return tdev(data, m; ci=ci, confidence=c)
        name === :hdev     && return hdev(data, m; ci=ci, confidence=c)
        name === :mhdev    && return mhdev(data, m; ci=ci, confidence=c)
        name === :htdev    && return htdev(data, m; ci=ci, confidence=c)
        name === :totdev   && return totdev(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :mtotdev  && return mtotdev(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :ttotdev  && return ttotdev(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :htotdev  && return htotdev(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :mhtotdev && return mhtotdev(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :theo1    && return theo1(data, m; ci=ci, confidence=c, correct_bias=b)
        name === :theoh    && return theoh(data, m; ci=ci, confidence=c)
        name === :pdev     && return pdev(data, m; ci=ci, confidence=c)
        name === :mtie     && return mtie(data, m; ci=ci, confidence=c)
        name === :tierms   && return tierms(data, m; ci=ci, confidence=c)
        throw(ArgumentError("unknown deviation " * String(name)))
    end
end

@static if !isdefined(SigmaTau, :tau_mode)
    function tau_mode(s::AbstractString)
        k = lowercase(s)
        k == "octave" && return Octave
        (k == "half-octave" || k == "halfoctave") && return HalfOctave
        (k == "quarter-octave" || k == "quarteroctave") && return QuarterOctave
        k == "decade" && return Decade
        (k == "half-decade" || k == "halfdecade") && return HalfDecade
        (k == "all" || k == "alltaus") && return AllTaus
        throw(ArgumentError("unknown τ grid \"" * String(s) * "\""))
    end
end

@static if !isdefined(SigmaTau, :summarize)
    # Gaps, 5σ outliers (indices) and the linear frequency drift per day.
    function summarize(data::PhaseData)
        x = data.x
        N = length(x)
        sxy = 0.0; sxx = 0.0; st = 0.0; sy = 0.0; n = 0
        for i in 1:N-1
            yi = (Float64(x[i + 1]) - Float64(x[i])) / data.tau0
            isfinite(yi) || continue
            ti = (i - 1) * data.tau0
            st += ti; sy += yi; sxy += ti * yi; sxx += ti * ti; n += 1
        end
        drift = n > 2 ? (n * sxy - st * sy) / (n * sxx - st * st) * 86400.0 : NaN
        return (n=N, tau0=Float64(data.tau0), span=N * data.tau0, gaps=count(isnan, x),
                outliers=find_outliers(data; nsigma=5), drift_per_day=drift)
    end
end

# The library's Welch segment once the registry lands is the largest power of
# two ≤ N/4 (at least 16); before that its default is 256, so pass it here.
function welch(f::F, data::PhaseData, nperseg::Int) where {F}
    nperseg > 0 && return f(data; nperseg=min(nperseg, length(data.x)))
    @static if isdefined(SigmaTau, :estimators)
        return f(data)
    else
        return f(data; nperseg=pow2_at_most(length(data.x) ÷ 4))
    end
end
