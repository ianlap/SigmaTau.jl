# Holdover test on a phase record: learn frequency (and drift) on a training
# stretch, free-run, and compare with what the clock did, at every start the
# record allows (Delporte et al. 2000).

"""
    HoldoverResult

Output of [`holdover`](@ref): the statistics of the time error at each
holdover time, over every start the record allows.

$(TYPEDFIELDS)
"""
struct HoldoverResult
    "Holdover times in seconds, from the last training sample."
    tau::Vector{Float64}
    "RMS time error (s)."
    rms::Vector{Float64}
    "The |error| (s) that `confidence` of the trials stay within."
    quantile::Vector{Float64}
    "Largest |error| (s) of any trial."
    worst::Vector{Float64}
    "Mean error (s); away from zero when the fit lacks a term the clock has, such as drift."
    mean::Vector{Float64}
    "Number of trials."
    n::Vector{Int}
    "Length of the training stretch in seconds."
    train::Float64
    "Whether the fit is a parabola (frequency drift) rather than a line."
    drift::Bool
    "Whether each trial starts from the last reading rather than the fit's value there."
    anchor::Bool
    "Level of `quantile`."
    confidence::Float64
    "The limit asked for (s); `NaN` when none."
    limit::Float64
    "Holdover time (s) at which `quantile` reaches `limit`; `NaN` when none."
    time_to_limit::Float64
end

Base.show(io::IO, h::HoldoverResult) =
    print(io, "HoldoverResult(", length(h.tau), " pts, τ∈[", _show4g(h.tau[1]), ", ",
          _show4g(h.tau[end]), "] s, trained on ", _show4g(h.train), " s)")

function Base.show(io::IO, ::MIME"text/plain", h::HoldoverResult)
    show(io, h)
    isnan(h.limit) || print(io, "\n  ", _show4g(h.limit), " s reached after ", _show4g(h.time_to_limit), " s")
    @printf(io, "\n  %10s  %10s  %10s  %10s  %8s", "tau", "rms", "p" * _show4g(100h.confidence), "worst", "n")
    for i in eachindex(h.tau)
        @printf(io, "\n  %10.4g  %10.4g  %10.4g  %10.4g  %8d", h.tau[i], h.rms[i], h.quantile[i], h.worst[i], h.n[i])
    end
end

# Tⱼ = Σᵢ (−i)ʲ v[s−i] over the n samples ending at s, j = 0 … length(T)−1:
# the moments of a window about its last sample.
function _window_moments!(T::Vector{Float64}, v::Vector{Float64}, s::Int, n::Int)
    fill!(T, 0.0)
    for i in 0:n-1, j in eachindex(T)
        T[j] += (-float(i))^(j - 1) * v[s-i]
    end
    return T
end

# The same moments after the window moves on one sample to end at s.
function _slide_moments!(T::Vector{Float64}, v::Vector{Float64}, s::Int, n::Int)
    for j in length(T)-1:-1:0
        t = 0.0
        for l in 0:j
            t += binomial(j, l) * (-1.0)^(j - l) * T[l+1]
        end
        T[j+1] = t - (-float(n))^j * v[s-n]
    end
    T[1] += v[s]
    return T
end

# Least-squares polynomial of degree d through each stretch of n samples of
# x, gaps skipped: row k holds the coefficients for the stretch ending at
# s = n + k − 1, in powers of the sample count from s. Rows are NaN where
# gaps take more than half the stretch.
function _sliding_fits(x::Vector{Float64}, n::Int, d::Int)
    N = length(x)
    w = Float64.(isfinite.(x))
    v = ifelse.(isfinite.(x), x, 0.0)
    W, Z = zeros(2d + 1), zeros(d + 1)
    scale = float(n) .^ (0:d)                       # the normal equations in units of the stretch
    A, C = zeros(d + 1, d + 1), fill(NaN, N - n, d + 1)
    for s in n:N-1
        # Summed afresh once per stretch length, so round-off cannot build up.
        if (s - n) % n == 0
            _window_moments!(W, w, s, n); _window_moments!(Z, v, s, n)
        else
            _slide_moments!(W, w, s, n); _slide_moments!(Z, v, s, n)
        end
        2W[1] >= n && W[1] > d || continue
        for j in 0:d, l in 0:d
            A[j+1, l+1] = W[j+l+1] / (scale[j+1] * scale[l+1])
        end
        C[s-n+1, :] = (_sym_inv(A) * (Z ./ scale)) ./ scale
    end
    return C
end

# Where a rising curve q(τ) reaches `limit`, interpolated on log axes: 0 when
# it starts above, Inf when it never gets there.
function _time_to_limit(τ::Vector{Float64}, q::Vector{Float64}, limit::Real)
    i = findfirst(>=(limit), q)
    i === nothing && return Inf
    i == 1 && return 0.0
    q[i-1] > 0 || return τ[i]
    return τ[i-1] * (τ[i] / τ[i-1])^(log(limit / q[i-1]) / log(q[i] / q[i-1]))
end

_whole_samples(t::Real, τ₀::Float64, what::String) = (k = round(Int, t / τ₀);
    k >= 1 && isapprox(k * τ₀, t; rtol=1e-9) ? k :
    throw(ArgumentError("holdover: $what = $t is not a whole multiple of tau0 = $τ₀")))

"""
    holdover(data::PhaseData; train, tau, drift=false, anchor=true, limit=nothing, confidence=0.95)
        → HoldoverResult

How far the clock wanders once its reference is gone, measured on the
record itself. A line is fitted to `train` seconds of phase (a parabola
with `drift = true`) and the clock is then left to run on the fit's
frequency, and drift, from the last reading of the stretch:

    e(τ) = x(t + τ) − x(t) − ŷ τ − ½ D̂ τ²

the time interval error of ITU-T G.810 less what the stretch taught. The
stretch then moves on one sample, and so on through the record, as in
Delporte et al. (2000). The result holds the RMS, the `confidence`
quantile, the worst case and the mean of e at each holdover time.

`train` defaults to a quarter of the record and `tau` to the octave
multiples of τ₀ up to `train`; both are whole multiples of τ₀. With `limit`
(s), `time_to_limit` is the holdover time at which the quantile reaches it:
`0` when it is passed at the first `tau` already, `Inf` when it is not
reached at the last.

With `anchor = false` the fit's own value at the last sample replaces the
reading x(t): the extrapolated fit of Delporte et al., which averages
white PM down but on a clock with random-walk FM starts from the fit's
misfit.

For a clock with white FM `q1` and random-walk FM `q2` (see
[`ClockModel`](@ref)), a line over T seconds gives the variance

    q1 τ (1 + 6τ/5T) + q2 (τ³/3 + 13 T τ²/35)

the first term of each bracket the clock's own noise over τ, the second
the error of the learned frequency. A line suits a clock without drift,
such as caesium, and then holds better than a parabola; with drift, the
line's mean error grows as τ².

The error includes the white PM of the readings. `NaN` samples are gaps: a
start is skipped when gaps take more than half its training stretch, the
reading it would be compared with or, when anchored, its last reading.

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 65536, 1.0; h = Dict(0 => 2e-22, -2 => 1e-30))
h = holdover(p; train = 4096, limit = 1e-9)
h.rms                 # at h.tau
h.time_to_limit       # seconds within 1 ns, in 95 % of the trials
```
"""
function holdover(data::PhaseData; train::Union{Real,Nothing}=nothing,
                  tau::Union{AbstractVector{<:Real},Nothing}=nothing, drift::Bool=false,
                  anchor::Bool=true, limit::Union{Real,Nothing}=nothing, confidence::Real=0.95)
    0 < confidence < 1 || throw(ArgumentError("holdover: confidence must be between 0 and 1, got $confidence"))
    τ₀, d = data.tau0, drift ? 2 : 1
    # The error of a degree-d fit is blind to a degree-d trend; taking the
    # record's own out first keeps the sums small.
    x = _detrend_core!(Float64.(data.x), drift ? :quadratic : :linear)
    N = length(x)
    n = train === nothing ? N ÷ 4 : _whole_samples(train, τ₀, "train")
    d + 1 < n < N || throw(ArgumentError(
        "holdover: a training stretch of $n samples does not fit a record of $N with $(d + 1) terms"))
    hs = tau === nothing ? _grid(Octave, min(n, (N - n) ÷ 2)) : [_whole_samples(t, τ₀, "tau") for t in tau]
    !isempty(hs) && issorted(hs) && hs[end] <= N - n || throw(ArgumentError(
        "holdover: tau must rise and stay within the $((N - n) * τ₀) s the record has past the training stretch"))
    C = _sliding_fits(x, n, d)
    rms, q, worst, μ = zeros(length(hs)), zeros(length(hs)), zeros(length(hs)), zeros(length(hs))
    cnt, e = zeros(Int, length(hs)), Float64[]
    for (i, h) in enumerate(hs)
        empty!(e)
        for s in n:N-h
            c = view(C, s - n + 1, :)
            y = x[s+h] - (anchor ? x[s] - c[1] : 0.0) - evalpoly(float(h), c)
            isfinite(y) && push!(e, y)
        end
        isempty(e) && throw(ArgumentError("holdover: gaps leave no trial at tau = $(h * τ₀)"))
        μ[i], rms[i], cnt[i] = sum(e) / length(e), sqrt(sum(abs2, e) / length(e)), length(e)
        e .= abs.(e)
        q[i], worst[i] = quantile(e, confidence), maximum(e)
    end
    τ = hs .* τ₀
    return HoldoverResult(τ, rms, q, worst, μ, cnt, n * τ₀, drift, anchor, Float64(confidence),
                          limit === nothing ? NaN : Float64(limit),
                          limit === nothing ? NaN : _time_to_limit(τ, q, limit))
end

# Weights on the phase samples 1 … n + h of the error `holdover` measures h
# samples past a stretch of n: +1 at n + h, less the prediction's.
function _holdover_weights(n::Int, h::Int, d::Int, anchor::Bool)
    u = collect(1-n:0) ./ n
    A = zeros(d + 1, d + 1)
    for j in 0:d, l in 0:d
        A[j+1, l+1] = sum(k -> k^(j + l), u)
    end
    Ainv, w = _sym_inv(A), zeros(n + h)
    for j in (anchor ? 1 : 0):d, l in 0:d
        w[1:n] .-= (h / n)^j * Ainv[j+1, l+1] .* u .^ l
    end
    anchor && (w[n] -= 1)
    w[n+h] += 1
    return w
end

# Σⱼₖ aⱼ aₖ / (1 − 4(j − k)²) = ∫₀^{π/2} sin θ |Σ aₖ e^{−2iθk}|² dθ, through the
# autocorrelation of a.
function _flicker_form(a::Vector{Float64})
    L = length(a)
    F = rfft(vcat(a, zeros(nextpow(2, 2L) - L)))
    r = irfft(abs2.(F), 2 * (length(F) - 1))
    return r[1] + 2 * sum(r[m+1] / (1 - 4.0 * m^2) for m in 1:L-1)
end

# Variance of the error with weights w under the unit component
# S_y(f) = (sin(πfτ₀)/(πτ₀))^α of a record sampled at τ₀ (as in
# `_powerlaw_var`). Summing w once gives the weights on the frequency
# samples, twice those on their steps, and so on; p sums leave the spectrum
# white (even α) or as sin θ (odd α), both exact in the time domain.
function _holdover_var(w::Vector{Float64}, α::Int, τ₀::Float64)
    p = (3 - α) ÷ 2
    a = w
    for _ in 1:p
        a = cumsum(a)
    end
    form = iseven(α) ? π / 2 * sum(abs2, a) : _flicker_form(a)
    return τ₀^2 * 4.0^(p - 1) / (π * τ₀)^(α + 1) * form
end

"""
    HoldoverPrediction

Output of [`holdover`](@ref) on a noise model: the time error the model
predicts at each holdover time, and what each noise type adds to it.

$(TYPEDFIELDS)
"""
struct HoldoverPrediction
    "Holdover times in seconds, from the last training sample."
    tau::Vector{Float64}
    "Standard deviation of the time error (s)."
    sigma::Vector{Float64}
    "The |error| (s) the clock stays within with probability `confidence`."
    quantile::Vector{Float64}
    "Exponents α of the model's noise types, highest first."
    alphas::Vector{Int}
    "Variance (s²) each noise type adds: one row per `tau`, one column per α. The rows sum to `sigma.^2`."
    budget::Matrix{Float64}
    "Length of the training stretch in seconds."
    train::Float64
    "Whether the fit is a parabola (frequency drift) rather than a line."
    drift::Bool
    "Whether the clock starts from the last reading rather than the fit's value there."
    anchor::Bool
    "Level of `quantile`."
    confidence::Float64
    "The limit asked for (s); `NaN` when none."
    limit::Float64
    "Holdover time (s) at which `quantile` reaches `limit`; `NaN` when none."
    time_to_limit::Float64
end

Base.show(io::IO, h::HoldoverPrediction) =
    print(io, "HoldoverPrediction(", length(h.tau), " pts, τ∈[", _show4g(h.tau[1]), ", ",
          _show4g(h.tau[end]), "] s, trained on ", _show4g(h.train), " s)")

function Base.show(io::IO, ::MIME"text/plain", h::HoldoverPrediction)
    show(io, h)
    isnan(h.limit) || print(io, "\n  ", _show4g(h.limit), " s reached after ", _show4g(h.time_to_limit), " s")
    @printf(io, "\n  %10s  %10s  %10s", "tau", "sigma", "p" * _show4g(100h.confidence))
    for α in h.alphas
        @printf(io, "  %10s", _nnls_symbol(α))
    end
    for i in eachindex(h.tau)
        @printf(io, "\n  %10.4g  %10.4g  %10.4g", h.tau[i], h.sigma[i], h.quantile[i])
        for j in eachindex(h.alphas)
            @printf(io, "  %10.4g", sqrt(h.budget[i, j]))
        end
    end
end

"""
    holdover(model; train, tau, tau0, drift=false, anchor=true, limit=nothing, confidence=0.95)
        → HoldoverPrediction

The same test predicted from a noise model: the standard deviation of the
error [`holdover`](@ref) measures on a record, for a clock with the
power-law noise of `model`. `model` is a [`PowerLawFit`](@ref) or a
[`ClockModel`](@ref) (whose `tau0` is the default), or a `Dict` of h_α
keyed by α or noise type, as in [`powerlaw_sigma`](@ref). `train` has no
default here.

The error is a weighted sum of phase samples, so its variance follows
exactly from the spectrum of the sampled record,
`h_α (sin(πfτ₀)/(πτ₀))^α` up to 1/(2τ₀), flicker noise included; `budget`
holds each noise type's share. For long stretches it gives the closed
forms of Delporte et al. (2000). The quantile is the Gaussian one,
1.96 σ at 95 %. A line takes α from 2 to −2; a parabola (`drift = true`)
also −3 and −4.

A drift D that the model does not hold adds the mean ½ D τ (τ + T) under a
line over T seconds, and nothing under a parabola.

# Examples

```julia
using SigmaTau
h = holdover(Dict(:WHFM => 2e-22, :FLFM => 1e-27); tau0 = 1.0, train = 3600, limit = 1e-9)
h.sigma               # at h.tau
sqrt.(h.budget)       # white FM and flicker FM, side by side
```
"""
function holdover(model; train::Real, tau::Union{AbstractVector{<:Real},Nothing}=nothing,
                  tau0::Real=_powerlaw_tau0(model), drift::Bool=false, anchor::Bool=true,
                  limit::Union{Real,Nothing}=nothing, confidence::Real=0.95)
    0 < confidence < 1 || throw(ArgumentError("holdover: confidence must be between 0 and 1, got $confidence"))
    τ₀, d = Float64(tau0), drift ? 2 : 1
    isfinite(τ₀) && τ₀ > 0 || throw(ArgumentError("holdover: pass `tau0`, the sample interval in seconds"))
    hd = _powerlaw_h(model)
    αs = reverse!(sort!(collect(keys(hd))))   # sort!(…; rev) is not trim-safe
    all(α -> -2d <= α <= 2, αs) || throw(ArgumentError(
        "holdover: α outside $(-2d)…2 has no finite error under a " * (drift ? "parabola" : "line; use drift = true")))
    n = _whole_samples(train, τ₀, "train")
    n > d + 1 || throw(ArgumentError("holdover: a training stretch of $n samples cannot hold $(d + 1) terms"))
    hs = tau === nothing ? _grid(Octave, n) : [_whole_samples(t, τ₀, "tau") for t in tau]
    !isempty(hs) && issorted(hs) || throw(ArgumentError("holdover: tau must rise"))
    budget = zeros(length(hs), length(αs))
    for (i, h) in enumerate(hs)
        w = _holdover_weights(n, h, d, anchor)
        for (j, α) in enumerate(αs)
            budget[i, j] = hd[α] * _holdover_var(w, α, τ₀)
        end
    end
    σ = sqrt.(vec(sum(budget; dims=2)))
    τ, q = hs .* τ₀, quantile(Normal(), (1 + confidence) / 2) .* σ
    return HoldoverPrediction(τ, σ, q, αs, budget, n * τ₀, drift, anchor, Float64(confidence),
                              limit === nothing ? NaN : Float64(limit),
                              limit === nothing ? NaN : _time_to_limit(τ, q, limit))
end
