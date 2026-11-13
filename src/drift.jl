# Frequency drift fits: a model is fitted by least squares to the fractional
# frequency y(t), and the drift rate is its slope at the end of the record
# (SP1065 §5.22; the Drift function of Stable32). A linear drift is also
# estimated by the methods of SP1065 §5.21 for phase data and by the
# four-point estimator of Greenhall (1997).

const DRIFT_MODELS = (:linear, :quadratic, :log, :diffusion)
const DRIFT_METHODS = (:frequency_fit, :phase_fit, :bisection, :second_difference, :four_point)

"""
    DriftFit

Result of [`drift_fit`](@ref): a frequency drift model fitted to ``y(t)``,
with the standard errors of its parameters and the drift rate at the end
of the record. The residuals are `f.y .- f.fit`, and `detrend(data, f)`
is the record without the fitted model.

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
    "Estimator of the parameters, one of the methods of [`drift_fit`](@ref)."
    method::Symbol
    "Parameters of the model, as in the table above."
    coef::Vector{Float64}
    "Standard error of each parameter; `NaN` where the method gives none."
    stderr::Vector{Float64}
    "Drift rate ``dy/dt`` at the end of the record, per day."
    rate::Float64
    "Standard error of `rate`, per day; `NaN` where the method gives none."
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
    drift_fit(data; model=:linear, method=:frequency_fit) → DriftFit

Fit of a frequency drift model to a record. The result gives the
parameters, their standard errors and the drift rate ``dy/dt`` at the end
of the record, per day.

# Models

- `:linear`, ``y(t) = a + bt``: a constant drift ``b``.
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

# Methods

`method` chooses the estimator of the linear model; the other models take
`:frequency_fit`. With ``N`` phase samples ``x_n`` at the interval ``τ_0``
and ``M = N - 1`` frequency samples ``y_n``:

| `method`             | Drift ``b``                                    | Least variance on |
|:---------------------|:-----------------------------------------------|:------------------|
| `:frequency_fit`     | slope of the least-squares line through ``y``  | white FM          |
| `:phase_fit`         | ``2c`` of the least-squares parabola ``a + bt + ct^2`` through ``x`` | white PM |
| `:bisection`         | difference of the mean frequencies of the two halves of the record, over half its length | |
| `:second_difference` | mean of the second differences of ``x``, over ``τ_0^2`` | random-walk FM |
| `:four_point`        | the four-point estimator, below                |                   |

These are the Linear Fit and Bisection Fit of Stable32 for frequency data
and its Quadratic Fit, Avg of 2nd Diffs and 3-Point Fit for phase data
([SP1065 §5.21 and §5.22](@cite riley-2008-sp1065)); the three-point fit
and the bisection are one estimator, ``4[x_N - 2x_{mid} + x_1]/T^2`` on a
record of length ``T`` with an odd number of phase samples
([Weiss and Hackman 1992, Eq. (1)](@cite weiss-1992-three-point)). The
four-point estimator of
[Greenhall 1997, Eq. (3)](@cite greenhall-1997-drift-estimator) is not in
Stable32. With ``w_n = \\sum_{j \\le n} x_j``, ``n_1`` the integer nearest
``N/10`` and ``r = n_1/N``, it is

```math
b = \\frac{6}{N^3 τ_0^2\\, r(1 - r)}
    \\left[w_N - w_0 - \\frac{w_{N-n_1} - w_{n_1}}{1 - 2r}\\right].
```

Every method returns the true drift on a record without noise. On a noisy
record the variance of the estimate depends on the noise, and
[Greenhall 1997, Table 1](@cite greenhall-1997-drift-estimator) gives it
for each method and power-law noise. The variance of the bisection is
within a factor 1.4 of the least on white FM, flicker FM and random-walk
FM. That of the four-point estimator is within a factor 1.6 of the least
on all five noises, white PM to random-walk FM. At a given record length,
the variances of the bisection and of the frequency fit on white PM, and
that of the mean second difference on white FM, grow with the number of
samples.

# Standard errors

The standard errors of `:frequency_fit` and `:phase_fit` are those of
least squares, from the scatter of the residuals, and that of
`:second_difference` is the standard deviation of the mean of the second
differences ([Barnes 1983, §II to §IV](@cite barnes-1983-drift)). Each
assumes the noise on which its method has the least variance. On a redder
noise it is too small, as on a record with flicker FM or random-walk FM:
on the quartz oscillators of Barnes 1983, Table 6, about 17 times for the
fit to phase and 5 times for the fit to frequency. On a whiter noise it is
too large. The deviation of the residuals shows whether they are white.

`:bisection` and `:four_point` have no standard error computed from the
record (`NaN`). On a record of length ``T`` their variances are

```math
\\frac{8}{T^2}\\,σ_y^2\\!\\left(\\frac{T}{2}\\right)
\\qquad\\text{and}\\qquad
\\frac{A^2}{T^2}\\,\\mathrm{Mod}\\,σ_y^2\\!\\left(\\frac{T}{3}\\right),
\\quad 3.14 \\le A \\le 3.80,
```

for any of the five noises
([Weiss and Hackman 1992, Eq. (7)](@cite weiss-1992-three-point);
[Greenhall 1997, §5](@cite greenhall-1997-drift-estimator)), with the
deviations of the record less its true drift. The record itself cannot
measure them at these averaging times, so they are extrapolated from the
deviations at shorter ones, by judgment of the noise type.

# Gaps

A `PhaseData` record is differenced to frequency first, except for
`:phase_fit` and `:four_point`, and missing samples are skipped.
`:four_point` needs a record without gaps, and so does `:phase_fit` of a
`FrequencyData` record, which is integrated to phase.

# Examples

A drift of 2 × 10⁻¹⁵ per second is 1.728 × 10⁻¹⁰ per day:

```jldoctest
julia> using SigmaTau

julia> y = FrequencyData(1e-9 .+ 2e-15 .* (0:999), 1.0);

julia> f = drift_fit(y);

julia> round(f.rate; sigdigits=4)
1.728e-10

julia> round(drift_fit(y; method=:four_point).rate; sigdigits=4)
1.728e-10
```
"""
function drift_fit(data::FrequencyData; model::Symbol=:linear, method::Symbol=:frequency_fit)
    return _drift_fit(Float64[], _f64(data.y), data.tau0, model, method)
end

function drift_fit(data::PhaseData; model::Symbol=:linear, method::Symbol=:frequency_fit)
    return _drift_fit(_f64(data.x), _f64(FrequencyData(data).y), data.tau0, model, method)
end

# `x` is the phase of a PhaseData record and empty for a FrequencyData one.
function _drift_fit(x::Vector{Float64}, y::Vector{Float64}, tau0::Float64, model::Symbol, method::Symbol)
    model in DRIFT_MODELS ||
        throw(ArgumentError("drift_fit: model must be one of $(join(DRIFT_MODELS, ", ")), got :$model"))
    method in DRIFT_METHODS ||
        throw(ArgumentError("drift_fit: method must be one of $(join(DRIFT_METHODS, ", ")), got :$method"))
    method === :frequency_fit || model === :linear ||
        throw(ArgumentError("drift_fit: method :$method estimates a linear drift; model :$model takes :frequency_fit"))
    t = collect(0:length(y)-1) .* tau0
    ok = isfinite.(y)
    nok = count(ok)
    p = model === :linear ? 2 : 3
    nok > p || throw(ArgumentError("drift_fit: need more than $p measured samples, got $nok"))
    coef, se, rate_se =
        method === :frequency_fit     ? _drift_frequency_fit(model, t, y, ok) :
        method === :phase_fit         ? _drift_phase_fit(_drift_phase(x, y, tau0, method), tau0) :
        method === :bisection         ? _drift_bisection(t, y, ok) :
        method === :second_difference ? _drift_second_difference(t, y, ok, tau0) :
                                        _drift_four_point(_drift_phase(x, y, tau0, method), t, y, tau0)
    fit = _drift_eval(model, coef, t)
    rss = sum(abs2, (y[i] - fit[i] for i in eachindex(y) if ok[i]))
    return DriftFit(model, method, coef, se, _drift_rate(model, coef, t[end]) * 86400.0, rate_se * 86400.0,
                    t, y, fit, sqrt(rss / nok))
end

# Least squares on the frequency: parameters, their standard errors and that
# of the drift rate at the end of the record (delta method).
function _drift_frequency_fit(model::Symbol, t::Vector{Float64}, y::Vector{Float64}, ok::BitVector)
    p = model === :linear ? 2 : 3
    coef = model === :log ? _drift_log(t, y, ok) :
           model === :diffusion ? _drift_diffusion(t, y, ok) : _drift_poly(t, y, ok, p)
    (model === :log || model === :diffusion) && (coef = _gauss_newton(model, coef, t, y, ok))
    cov = _drift_cov(model, coef, t, y, ok)
    g = _drift_rate_grad(model, coef, t[end])
    rate_var = sum(g[i] * cov[i, j] * g[j] for i in 1:p, j in 1:p)
    return coef, [sqrt(max(cov[i, i], 0.0)) for i in 1:p], sqrt(max(rate_var, 0.0))
end

# Covariance of the least-squares parameters, from the scatter of the residuals.
function _drift_cov(model::Symbol, coef::Vector{Float64}, t::Vector{Float64}, y::Vector{Float64}, ok::BitVector)
    fit = _drift_eval(model, coef, t)
    rss = sum(abs2, (y[i] - fit[i] for i in eachindex(y) if ok[i]))
    return (rss / (count(ok) - length(coef))) .* _inv_normal(_drift_jac(model, coef, t), ok)
end

# Phase for the methods that work on it: the record's own, or the frequency
# integrated, which a gap prevents.
function _drift_phase(x::Vector{Float64}, y::Vector{Float64}, tau0::Float64, method::Symbol)
    isempty(x) || return x
    all(isfinite, y) ||
        throw(ArgumentError("drift_fit: method :$method on a frequency record needs one without gaps; fill them with fillgaps first"))
    return _freq_to_phase(FrequencyData(y, tau0)).x
end

# Least-squares parabola a + bt + ct² through the phase (optimum on white PM;
# Barnes 1983, App. A), as the line it gives the frequency samples: these
# are the first differences of the parabola, b + cτ₀ + 2ct.
function _drift_phase_fit(x::Vector{Float64}, tau0::Float64)
    t = collect(0:length(x)-1) .* tau0
    ok = isfinite.(x)
    nok = count(ok)
    nok > 3 || throw(ArgumentError("drift_fit: method :phase_fit needs more than 3 measured phase samples, got $nok"))
    c = _drift_poly(t, x, ok, 3)
    cov = _drift_cov(:quadratic, c, t, x, ok)
    var_a = cov[2, 2] + 2tau0 * cov[2, 3] + tau0^2 * cov[3, 3]
    se = 2sqrt(max(cov[3, 3], 0.0))
    return [c[2] + c[3] * tau0, 2c[3]], [sqrt(max(var_a, 0.0)), se], se
end

# Mean of the measured samples of y over the index range r, and their mean time.
function _range_mean(t::Vector{Float64}, y::Vector{Float64}, ok::BitVector, r::UnitRange{Int})
    n = 0
    sy = st = 0.0
    @inbounds for i in r
        ok[i] || continue
        n += 1; sy += y[i]; st += t[i]
    end
    return sy / n, st / n
end

# Difference of the mean frequencies of the two halves of the record over
# the distance between their mean times: the line through the two means
# (SP1065 §5.22; on phase, the three-point estimator of Weiss and Hackman
# 1992, Eq. (1)). The middle sample of an odd count is left out.
function _drift_bisection(t::Vector{Float64}, y::Vector{Float64}, ok::BitVector)
    M = length(y)
    h = M ÷ 2
    y1, t1 = _range_mean(t, y, ok, 1:h)
    y2, t2 = _range_mean(t, y, ok, M-h+1:M)
    isfinite(y1) && isfinite(y2) ||
        throw(ArgumentError("drift_fit: method :bisection needs measured samples in both halves of the record"))
    b = (y2 - y1) / (t2 - t1)
    return [y1 - b * t1, b], [NaN, NaN], NaN
end

# Mean of the first differences of the frequency (the second differences of
# the phase) over τ₀, with the standard deviation of the mean (optimum on
# random-walk FM; Barnes 1983, §IV). The line starts at the first measured
# sample; without gaps it ends at the last.
function _drift_second_difference(t::Vector{Float64}, y::Vector{Float64}, ok::BitVector, tau0::Float64)
    n = 0
    d = 0.0
    @inbounds for i in 1:length(y)-1
        ok[i] && ok[i+1] || continue
        n += 1; d += y[i+1] - y[i]
    end
    n > 0 || throw(ArgumentError("drift_fit: method :second_difference needs two measured samples in a row"))
    d /= n
    ss = 0.0
    @inbounds for i in 1:length(y)-1
        ok[i] && ok[i+1] && (ss += abs2(y[i+1] - y[i] - d))
    end
    se = n > 1 ? sqrt(ss / ((n - 1) * n)) / tau0 : NaN
    b = d / tau0
    i1 = findfirst(ok)::Int
    return [y[i1] - b * t[i1], b], [NaN, se], se
end

# Four-point estimator of Greenhall (1997), Eq. (3), from the cumulative sum
# w of the phase at 0, n₁, N − n₁ and N. w(N) − w(0) is the sum of the phase
# and w(N − n₁) − w(n₁) that of its middle. The estimator ignores a line in
# the phase, so the one through the end points is taken out first, which
# keeps the sums small. The frequency line has the mean of the record.
function _drift_four_point(x::Vector{Float64}, t::Vector{Float64}, y::Vector{Float64}, tau0::Float64)
    all(isfinite, x) || throw(ArgumentError("drift_fit: method :four_point needs a record without gaps; fill them with fillgaps first"))
    N = length(x)
    n1 = max(1, round(Int, N / 10))
    r = n1 / N
    z = _detrend_core!(copy(x), :endpoint)
    b = 6 / (float(N)^3 * tau0^2 * r * (1 - r)) * (sum(z) - sum(view(z, n1+1:N-n1)) / (1 - 2r))
    return [sum(y) / length(y) - b * sum(t) / length(t), b], [NaN, NaN], NaN
end

"""
    detrend(data, fit::DriftFit) → PhaseData or FrequencyData

Copy of a record with the frequency model of its [`drift_fit`](@ref)
removed: `fit.fit` from the frequency samples, or its running sum times
``τ_0`` from the phase samples, which leaves the first one as it is.

The deviation of the result shows whether the model left white residuals,
and it is the deviation that the variance of the drift estimate is read
from ([Weiss and Hackman 1992](@cite weiss-1992-three-point)). Removing an
estimated drift lowers the deviations at the longest averaging times:
after `:bisection` on an odd number of phase samples, the Allan deviation
at half the record is zero.

# Examples

```jldoctest
julia> using SigmaTau

julia> x = PhaseData([2.0 + 3k + 0.5k^2 for k in 0:10], 1.0);

julia> r = detrend(x, drift_fit(x; method=:bisection));

julia> maximum(abs, r.x .- 2.0) < 1e-12
true
```
"""
function detrend(fd::FrequencyData, f::DriftFit)
    _check_fit(f, length(fd.y), fd.tau0)
    return FrequencyData(_f64(fd.y) .- f.fit, fd.tau0; source=fd.source, start=fd.start)
end

function detrend(pd::PhaseData, f::DriftFit)
    _check_fit(f, length(pd.x) - 1, pd.tau0)
    x = Vector{Float64}(pd.x)
    acc = 0.0
    @inbounds for i in eachindex(f.fit)
        acc += f.fit[i]
        x[i+1] -= acc * pd.tau0
    end
    return PhaseData(x, pd.tau0; source=pd.source, start=pd.start)
end

function _check_fit(f::DriftFit, M::Int, tau0::Float64)
    length(f.fit) == M && (M < 2 || f.t[2] == tau0) ||
        throw(ArgumentError("detrend: the fit is of another record (its length or sample interval differs)"))
    return nothing
end

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

# " ± s", or nothing where the method gives no standard error.
_show_pm(s::Float64) = isnan(s) ? "" : " ± " * _show4g(s)

function Base.show(io::IO, f::DriftFit)
    print(io, "DriftFit(:", f.model, ", :", f.method, ", rate=", _show4g(f.rate), _show_pm(f.rate_stderr), " /day)")
end

function Base.show(io::IO, ::MIME"text/plain", f::DriftFit)
    names = f.model === :linear ? ("a", "b") : ("a", "b", "c")
    print(io, "DriftFit: ", f.model, " model by ", f.method, ", ", count(isfinite, f.y),
          " samples, rms residual ", _show4g(f.rms))
    for (n, c, s) in zip(names, f.coef, f.stderr)
        print(io, "\n  ", n, " = ", _show4g(c), _show_pm(s))
    end
    print(io, "\n  drift at end = ", _show4g(f.rate), _show_pm(f.rate_stderr), " per day")
end
