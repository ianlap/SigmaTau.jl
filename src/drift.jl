# Frequency drift fits: a model is fitted by least squares to the fractional
# frequency y(t), and the drift rate is its slope at the end of the record
# (SP1065 §5.22; the Drift function of Stable32).

const DRIFT_MODELS = (:linear, :quadratic, :log, :diffusion)

"""
    DriftFit

Result of [`drift_fit`](@ref): a frequency drift model fitted to ``y(t)``,
with the standard errors of its parameters and the drift rate at the end
of the record. The residuals are `f.y .- f.fit`.

`coef` holds the parameters of the model, with ``t`` in seconds from the
first sample:

| `model`      | ``y(t)``               | `coef`      |
|:-------------|:-----------------------|:------------|
| `:linear`    | ``a + bt``             | `[a, b]`    |
| `:quadratic` | ``a + bt + ct^2``      | `[a, b, c]` |
| `:log`       | ``a\\ln(1 + bt) + c``   | `[a, b, c]` |
| `:diffusion` | ``a + b\\sqrt{t + c}``  | `[a, b, c]` |

# Fields

$(TYPEDFIELDS)
"""
struct DriftFit
    "Model fitted: `:linear`, `:quadratic`, `:log` or `:diffusion`."
    model::Symbol
    "Parameters of the model, as in the table above."
    coef::Vector{Float64}
    "Standard error of each parameter."
    stderr::Vector{Float64}
    "Drift rate ``dy/dt`` at the end of the record, per day."
    rate::Float64
    "Standard error of `rate`, per day."
    rate_stderr::Float64
    "Sample times, in seconds from the first sample."
    t::Vector{Float64}
    "Fractional frequency; `NaN` at a missing sample."
    y::Vector{Float64}
    "The model at each `t`."
    fit::Vector{Float64}
    "Rms of the residuals `y − fit` over the measured samples."
    rms::Float64
end

"""
    drift_fit(data; model=:linear) → DriftFit

Least-squares fit of a frequency drift model to the fractional frequency
of a record. The result gives the parameters, their standard errors and
the drift rate ``dy/dt`` at the end of the record, per day.

# Models

- `:linear`, ``y(t) = a + bt``: a constant drift ``b``, the optimum fit
  under white FM noise.
- `:quadratic`, ``y(t) = a + bt + ct^2``: a drift that changes linearly
  with time.
- `:log`, ``y(t) = a\\ln(1 + bt) + c``: the stabilization of an oscillator.
- `:diffusion`, ``y(t) = a + b\\sqrt{t + c}``.

The linear, log and diffusion models are the frequency drift models of
[SP1065 §5.22](@cite riley-2008-sp1065) and of the Drift function of
[Stable32](@cite hamilton-2008-stable32-manual); the constant ``c`` of the
log model is in the form its test suite reports. ``t`` runs from 0 at the
first sample, so ``a`` of the linear model is the frequency there.
Stable32 counts the samples from 1 and reports the intercept one sample
interval earlier; the slope is the same.

A `PhaseData` record is differenced to frequency first, and missing
samples are skipped.

The standard errors are those of least squares, from the scatter of the
residuals. They assume white residuals and are too small where the
residuals are not white, as on a record with flicker FM or random-walk FM.
The deviation of the residuals shows whether they are.

# Examples

A drift of 2 × 10⁻¹⁵ per second is 1.728 × 10⁻¹⁰ per day:

```jldoctest
julia> using SigmaTau

julia> y = FrequencyData(1e-9 .+ 2e-15 .* (0:999), 1.0);

julia> f = drift_fit(y);

julia> round(f.rate; sigdigits=4)
1.728e-10
```
"""
function drift_fit(data::FrequencyData; model::Symbol=:linear)
    model in DRIFT_MODELS ||
        throw(ArgumentError("drift_fit: model must be one of $(join(DRIFT_MODELS, ", ")), got :$model"))
    y = _f64(data.y)
    t = collect(0:length(y)-1) .* data.tau0
    ok = isfinite.(y)
    nok = count(ok)
    p = model === :linear ? 2 : 3
    nok > p || throw(ArgumentError("drift_fit: need more than $p measured samples, got $nok"))
    coef = model === :log ? _drift_log(t, y, ok) :
           model === :diffusion ? _drift_diffusion(t, y, ok) : _drift_poly(t, y, ok, p)
    (model === :log || model === :diffusion) && (coef = _gauss_newton(model, coef, t, y, ok))
    J = _drift_jac(model, coef, t)
    fit = _drift_eval(model, coef, t)
    rss = sum(abs2, (y[i] - fit[i] for i in eachindex(y) if ok[i]))
    s2 = rss / (nok - p)
    cov = s2 .* _inv_normal(J, ok)
    g = _drift_rate_grad(model, coef, t[end])
    rate_var = sum(g[i] * cov[i, j] * g[j] for i in 1:p, j in 1:p)
    return DriftFit(model, coef, [sqrt(max(cov[i, i], 0.0)) for i in 1:p],
                    _drift_rate(model, coef, t[end]) * 86400.0, sqrt(max(rate_var, 0.0)) * 86400.0,
                    t, y, fit, sqrt(rss / nok))
end

drift_fit(data::PhaseData; kwargs...) = drift_fit(FrequencyData(data); kwargs...)

function _drift_eval(model::Symbol, c::Vector{Float64}, t::Vector{Float64})
    model === :linear    && return c[1] .+ c[2] .* t
    model === :quadratic && return c[1] .+ c[2] .* t .+ c[3] .* t .^ 2
    model === :log       && return c[1] .* log1p.(c[2] .* t) .+ c[3]
    return c[1] .+ c[2] .* sqrt.(t .+ c[3])
end

function _drift_rate(model::Symbol, c::Vector{Float64}, t::Float64)
    model === :linear    && return c[2]
    model === :quadratic && return c[2] + 2c[3] * t
    model === :log       && return c[1] * c[2] / (1 + c[2] * t)
    return c[2] / (2sqrt(t + c[3]))
end

# ∂y/∂coef at each t.
function _drift_jac(model::Symbol, c::Vector{Float64}, t::Vector{Float64})
    one_ = ones(length(t))
    model === :linear    && return hcat(one_, t)
    model === :quadratic && return hcat(one_, t, t .^ 2)
    model === :log       && return hcat(log1p.(c[2] .* t), c[1] .* t ./ (1 .+ c[2] .* t), one_)
    return hcat(one_, sqrt.(t .+ c[3]), c[2] ./ (2 .* sqrt.(t .+ c[3])))
end

# Polish a nonlinear fit with Gauss–Newton steps on all three parameters;
# the profile search leaves the flat direction good to only ~1e-7.
function _gauss_newton(model::Symbol, c::Vector{Float64}, t, y, ok)
    rss(c) = sum(abs2(y[i] - f) for (i, f) in enumerate(_drift_eval(model, c, t)) if ok[i])
    best = rss(c)
    for _ in 1:8
        r = y .- _drift_eval(model, c, t)
        step, _ = _lsq(_drift_jac(model, c, t), r, ok)
        trial = c .+ step
        (model === :log ? trial[2] > 0 : trial[3] > -t[1]) || break
        v = rss(trial)
        # Near the optimum the residual is flat to rounding, so a step that
        # leaves it unchanged to 1e-12 still sharpens the parameters.
        v <= best * (1 + 1e-12) || break
        c, best = trial, min(v, best)
        maximum(abs.(step) ./ max.(abs.(c), floatmin())) < 1e-13 && break
    end
    return c
end

# ∂rate/∂coef at time t, for the delta-method standard error.
function _drift_rate_grad(model::Symbol, c::Vector{Float64}, t::Float64)
    model === :linear    && return [0.0, 1.0]
    model === :quadratic && return [0.0, 1.0, 2t]
    if model === :log
        u = 1 + c[2] * t
        return [c[2] / u, c[1] / u^2, 0.0]
    end
    s = sqrt(t + c[3])
    return [0.0, 1 / (2s), -c[2] / (4s^3)]
end

# (JᵀJ)⁻¹ over the measured rows, for the parameter covariance; the fit
# itself is `_lsq` (src/io/detrend.jl).
function _inv_normal(J::Matrix{Float64}, ok::BitVector)
    G = _gram(J, ok)
    sc = _col_scale(G)
    Ainv = _normal_inv(G, sc)
    return [Ainv[j, k] / (sc[j] * sc[k]) for j in axes(Ainv, 1), k in axes(Ainv, 2)]
end

# Polynomials are fitted in t/T (T the record span) for conditioning and
# converted back to seconds.
function _drift_poly(t, y, ok, p)
    T = max(t[end], 1.0)
    u = t ./ T
    X = p == 2 ? hcat(ones(length(t)), u) : hcat(ones(length(t)), u, u .^ 2)
    beta, _ = _lsq(X, y, ok)
    return p == 2 ? [beta[1], beta[2] / T] : [beta[1], beta[2] / T, beta[3] / T^2]
end

# The nonlinear models are linear in two parameters once the third is
# fixed: minimize the residual over that one (b for :log, c for :diffusion)
# on a log grid, then refine by golden section.
function _profile_min(rss::Function, lo::Float64, hi::Float64)
    grid = range(lo, hi; length=81)
    vals = [rss(g) for g in grid]
    k = argmin(vals)
    a, b = grid[max(k - 1, 1)], grid[min(k + 1, length(grid))]
    phi = (sqrt(5) - 1) / 2
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = rss(c), rss(d)
    for _ in 1:80
        if fc < fd
            b, d, fd = d, c, fc
            c = b - phi * (b - a); fc = rss(c)
        else
            a, c, fc = c, d, fd
            d = a + phi * (b - a); fd = rss(d)
        end
    end
    return (a + b) / 2
end

function _drift_log(t, y, ok)
    T = max(t[end], 1.0)
    one_ = ones(length(t))
    rss(lb) = last(_lsq(hcat(log1p.((10.0^lb / T) .* t), one_), y, ok))
    b = 10.0^_profile_min(rss, -3.0, 4.0) / T
    a, c = first(_lsq(hcat(log1p.(b .* t), one_), y, ok))
    return [a, b, c]
end

function _drift_diffusion(t, y, ok)
    T = max(t[end], 1.0)
    one_ = ones(length(t))
    rss(lc) = last(_lsq(hcat(one_, sqrt.(t .+ 10.0^lc * T)), y, ok))
    c = 10.0^_profile_min(rss, -4.0, 3.0) * T
    a, b = first(_lsq(hcat(one_, sqrt.(t .+ c)), y, ok))
    return [a, b, c]
end

function Base.show(io::IO, f::DriftFit)
    print(io, "DriftFit(:", f.model, ", rate=", _show4g(f.rate), " ± ", _show4g(f.rate_stderr), " /day)")
end

function Base.show(io::IO, ::MIME"text/plain", f::DriftFit)
    names = f.model === :linear ? ("a", "b") : ("a", "b", "c")
    print(io, "DriftFit: ", f.model, " model, ", count(isfinite, f.y), " samples, rms residual ",
          _show4g(f.rms))
    for (n, c, s) in zip(names, f.coef, f.stderr)
        print(io, "\n  ", n, " = ", _show4g(c), " ± ", _show4g(s))
    end
    print(io, "\n  drift at end = ", _show4g(f.rate), " ± ", _show4g(f.rate_stderr), " per day")
end
