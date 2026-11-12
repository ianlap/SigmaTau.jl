# Single-clock Kalman filter on a phase record: filter, Rauch–Tung–Striebel
# smoother, holdover error and an innovation check, on a ClockModel.

"""
    KalmanResult

Output of [`kalman_filter`](@ref) or [`kalman_smooth`](@ref): the state
estimate at every sample, its covariance, and the innovations.

# Fields

$(TYPEDFIELDS)
"""
struct KalmanResult
    "Clock model the filter ran on."
    model::ClockModel
    "Sample times in seconds from the first sample."
    t::Vector{Float64}
    "State estimates, one row per sample: phase (s), frequency, drift (1/s) for a model with a drift state, then the Markov frequency states, which add to the frequency, and the Markov phase states, which add to the readings."
    x::Matrix{Float64}
    "State covariances, `P[:, :, k]` at sample k."
    P::Array{Float64,3}
    "Innovations ``ν_k``, each reading less its prediction (s); `NaN` at gaps."
    innovation::Vector{Float64}
    "Innovation variances ``S_k`` (s²); `NaN` at gaps."
    innovation_var::Vector{Float64}
    "Whether `x` and `P` are smoothed."
    smoothed::Bool
    # Predicted states and covariances, kept for the smoother.
    xpred::Matrix{Float64}
    Ppred::Array{Float64,3}
end

function Base.show(io::IO, ::MIME"text/plain", r::KalmanResult)
    print(io, r.smoothed ? "Smoothed " : "", "KalmanResult, ", length(r.t), " samples, ",
          size(r.x, 2), " states")
end

# A / B for a covariance B whose states differ in scale by many decades
# (phase in s, drift in 1/s): solve in units of B's own standard deviations.
function _rdiv_cov(A::AbstractMatrix, B::AbstractMatrix)
    d = sqrt.(max.(diag(B), floatmin()))
    return ((A ./ d') / (B ./ (d * d'))) ./ d'
end

# Start from a least-squares polynomial through the first samples, so the
# filter needs no prior from the user.
function _kalman_start(z::Vector{Float64}, τ₀::Float64, n::Int)
    idx = findall(isfinite, z)
    k = idx[1:min(length(idx), max(4n, 16))]
    length(k) > n || throw(ArgumentError("kalman_filter: needs more than $n finite samples"))
    t = (k .- k[1]) .* τ₀
    # Filled in place: splatting a generator into hcat is not trim-safe.
    A = Matrix{Float64}(undef, length(k), n)
    for j in 1:n
        A[:, j] .= t .^ (j - 1) ./ factorial(j - 1)
    end
    c = A \ z[k]
    s² = max(sum(abs2, z[k] - A * c) / (length(k) - n), eps(maximum(abs, z[k]))^2)
    C = s² .* inv(A' * A)
    # Carry the estimate at t[1] back to sample 1 with zero added noise.
    Φ = [j >= i ? (-(k[1] - 1) * τ₀)^(j - i) / factorial(j - i) : 0.0 for i in 1:n, j in 1:n]
    return Φ * c, 10 .* (Φ * C * Φ')
end

# Stationary variance of each Markov state, frequency states first.
_markov_vars(m::ClockModel) = vcat(fill(_fm_var(m), length(_fm_taus(m))), fill(_pm_var(m), length(_pm_taus(m))))

# The start for model `m`: the polynomial start for phase, frequency and
# drift; the Markov states at zero with their stationary variance. The
# readings the phase starts from carry the Markov phase states. With
# `zero_drift` the drift starts at zero with the variance the parabola gives
# it, and phase and frequency start from a line: for a filter that cannot
# observe the drift itself, only its differences between clocks.
function _model_start(z::Vector{Float64}, m::ClockModel; widen::Float64=1.0, zero_drift::Bool=false)
    np, n = _npoly(m), _nstates(m)
    xp, Pp = _kalman_start(z, m.tau0, np)
    x, P = zeros(n), zeros(n, n)
    if zero_drift && np == 3
        x[1:2], P[1:2, 1:2] = _kalman_start(z, m.tau0, 2)
        P[3, 3] = Pp[3, 3]
    else
        x[1:np] = xp
        P[1:np, 1:np] = widen .* Pp
    end
    for (j, v) in enumerate(_markov_vars(m))
        P[np+j, np+j] = v
    end
    for j in n-length(_pm_taus(m))+1:n
        P[1, 1] += P[j, j]
        P[1, j] = P[j, 1] = -P[j, j]
    end
    return x, P
end

"""
    kalman_filter(data::PhaseData, m::ClockModel; x0=nothing, P0=nothing)
        → KalmanResult

Kalman filter of clock model `m` over a phase record.

Each sample is a reading ``z_k = h^{\\mathsf T} X_k + v_k`` of the state
``X_k`` of [`ClockModel`](@ref): the phase, plus the Markov phase states of
a model with flicker PM, plus white noise ``v_k`` of variance ``R``. From
one sample to the next the state moves by ``Φ``
([`state_transition`](@ref)) and gathers noise of covariance ``Q``
([`process_noise`](@ref)). The filter is the standard recursion
([Särkkä and Svensson (2023), Theorem 6.6](@cite sarkka-2023-bayesian-filtering)):

```math
\\begin{aligned}
X^-_k &= Φ X_{k-1}, & P^-_k &= Φ P_{k-1} Φ^{\\mathsf T} + Q, \\\\
ν_k &= z_k - h^{\\mathsf T} X^-_k, & S_k &= h^{\\mathsf T} P^-_k h + R, \\\\
K_k &= P^-_k h / S_k, & X_k &= X^-_k + K_k ν_k .
\\end{aligned}
```

The covariance is updated in the Joseph form
([Stein (2003), Eq. (5.8)](@cite stein-2003-demystified)), which stays
symmetric and positive semidefinite under round-off:

```math
P_k = (I - K_k h^{\\mathsf T})\\, P^-_k\\, (I - K_k h^{\\mathsf T})^{\\mathsf T}
      + R\\, K_k K_k^{\\mathsf T} .
```

A `NaN` sample is a gap: the state is predicted across it without an
update.

The filter needs no prior from the caller. It starts from a least-squares
line through the first 16 readings (a parabola for a model with drift),
with the covariance of that fit multiplied by ten, and with the Markov
states at zero with their stationary variances.

Pass the result to [`kalman_smooth`](@ref), [`holdover_error`](@ref),
[`innovation_check`](@ref) or [`prediction_error`](@ref).

# Arguments

- `data`: the phase record. Its `tau0` must be the model's.
- `m`: the clock model. A model with no noise at all is refused: the gain
  would be undefined once the start is used up.

# Keywords

- `x0`, `P0`: state and covariance at the first sample, before its reading
  is used, in place of the default start.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(2e-9 .* (0:199), 1.0);    # a clock that gains 2 ns each second

julia> m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26);

julia> kf = kalman_filter(p, m)
KalmanResult, 200 samples, 2 states

julia> kf.x[end, 2] ≈ 2e-9
true
```
"""
function kalman_filter(data::PhaseData, m::ClockModel;
                       x0::Union{AbstractVector,Nothing}=nothing,
                       P0::Union{AbstractMatrix,Nothing}=nothing)
    isapprox(data.tau0, m.tau0; rtol=1e-9) || throw(ArgumentError(
        "kalman_filter: record τ₀ = $(data.tau0) s but the model's is $(m.tau0) s"))
    _noiseless(m) && throw(ArgumentError("kalman_filter: the clock model has no noise, so the filter cannot weigh a reading"))
    z = Float64.(data.x)
    N, n = length(z), _nstates(m)
    xs, Ps0 = x0 === nothing || P0 === nothing ? _model_start(z, m) : (x0, P0)
    x̂ = Float64.(x0 === nothing ? xs : x0)
    P = Float64.(P0 === nothing ? Ps0 : P0)
    size(P) == (n, n) && length(x̂) == n || throw(ArgumentError("kalman_filter: x0 and P0 must have $n states"))
    Φ, Q = state_transition(m, m.tau0), process_noise(m, m.tau0)
    h, R = _reading_row(m), _reading_var(m)
    X, Xp = zeros(N, n), zeros(N, n)
    Ps, Pps = zeros(n, n, N), zeros(n, n, N)
    ν, S = fill(NaN, N), fill(NaN, N)
    for k in 1:N
        if k > 1
            x̂ = Φ * x̂
            P = Φ * P * Φ' + Q
        end
        Xp[k, :] = x̂; Pps[:, :, k] = P
        if isfinite(z[k])
            Ph = P * h
            S[k] = h' * Ph + R
            ν[k] = z[k] - h' * x̂
            K = Ph ./ S[k]
            x̂ = x̂ .+ K .* ν[k]
            IKH = Matrix{Float64}(I, n, n) .- K * h'
            P = IKH * P * IKH' .+ R .* (K * K')
            P = (P + P') ./ 2
        end
        X[k, :] = x̂; Ps[:, :, k] = P
    end
    t = (0:N-1) .* m.tau0
    return KalmanResult(m, collect(t), X, Ps, ν, S, false, Xp, Pps)
end

# The clock-model functions follow the phase; a frequency record is refused by name.
_needs_phase(f::AbstractString) = throw(ArgumentError(
    "$f: takes a phase record; convert the frequency record first with PhaseData(y)"))

kalman_filter(::FrequencyData, ::ClockModel; kw...) = _needs_phase("kalman_filter")

"""
    kalman_smooth(kf::KalmanResult) → KalmanResult

Rauch–Tung–Striebel smoother over a filtered run: every estimate then uses
the whole record, the samples after it as well as those before.

Running backward from the last sample
([Särkkä and Svensson (2023), Theorem 12.2](@cite sarkka-2023-bayesian-filtering)),

```math
\\begin{aligned}
G_k &= P_k Φ^{\\mathsf T} (P^-_{k+1})^{-1}, \\\\
X^s_k &= X_k + G_k (X^s_{k+1} - X^-_{k+1}), \\\\
P^s_k &= P_k + G_k (P^s_{k+1} - P^-_{k+1})\\, G_k^{\\mathsf T} ,
\\end{aligned}
```

with ``X_k``, ``P_k`` the filter's estimates and ``X^-_k``, ``P^-_k`` its
predictions. The innovations of the result are the filter's.

# Examples

```jldoctest
julia> using SigmaTau

julia> m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26);

julia> kf = kalman_filter(PhaseData(zeros(200), 1.0), m);

julia> ks = kalman_smooth(kf)
Smoothed KalmanResult, 200 samples, 2 states

julia> ks.P[1, 1, 100] < kf.P[1, 1, 100]    # the later samples add information
true
```
"""
function kalman_smooth(kf::KalmanResult)
    kf.smoothed && throw(ArgumentError("kalman_smooth: the result is already smoothed"))
    X, P = copy(kf.x), copy(kf.P)
    Φ = state_transition(kf.model, kf.model.tau0)
    for k in size(X, 1)-1:-1:1
        C = _rdiv_cov(kf.P[:, :, k] * Φ', kf.Ppred[:, :, k+1])
        X[k, :] = kf.x[k, :] .+ C * (X[k+1, :] .- kf.xpred[k+1, :])
        Pk = kf.P[:, :, k] .+ C * (P[:, :, k+1] .- kf.Ppred[:, :, k+1]) * C'
        P[:, :, k] = (Pk + Pk') ./ 2
    end
    return KalmanResult(kf.model, kf.t, X, P, kf.innovation, kf.innovation_var, true, kf.xpred, kf.Ppred)
end

"""
    holdover_error(kf::KalmanResult, tau) → Vector{Float64}

Standard deviation, in seconds, of the clock's time error after running
free for each τ in `tau` (seconds) from the last sample, as the model gives
it:

```math
σ_x(τ) = \\sqrt{\\left[Φ(τ)\\, P\\, Φ(τ)^{\\mathsf T} + Q(τ)\\right]_{11}} ,
```

with ``P`` the covariance of the filter at the last sample. The first term
is the error of the last estimate carried forward and the second the noise
the clock gathers on the way. It is the error of the clock's phase: the
noise of a reading, white or flicker PM, is not in it. [`holdover`](@ref)
measures the same error on the record.

# Examples

```jldoctest
julia> using SigmaTau

julia> m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26);

julia> kf = kalman_filter(PhaseData(zeros(200), 1.0), m);

julia> round.(holdover_error(kf, [10.0, 100.0, 1000.0]); sigdigits = 3)
3-element Vector{Float64}:
 4.82e-11
 1.63e-10
 2.12e-9
```
"""
function holdover_error(kf::KalmanResult, tau::AbstractVector{<:Real})
    P = kf.P[:, :, end]
    return map(tau) do τ
        Φ = state_transition(kf.model, τ)
        sqrt((Φ * P * Φ' + process_noise(kf.model, τ))[1, 1])
    end
end

"""
    innovation_check(kf::KalmanResult; confidence=0.95) → NamedTuple

Test of the model against the record, from the innovations of the filter.

Under the model the normalized innovations ``e_k = ν_k/\\sqrt{S_k}`` are
independent, with zero mean and unit variance
([Särkkä and Svensson (2023), Theorem 6.6](@cite sarkka-2023-bayesian-filtering)).
Two statistics of the ``n`` innovations test that:

- The sum of ``e_k^2`` is then χ² distributed with ``n`` degrees of
  freedom, which bounds their mean. A mean above the interval says the
  noise of the model is too low; below, too high.
- The lag-1 autocorrelation of ``e_k`` then stays within
  ``± z_{1-α/2}/\\sqrt{n}``, the band of an uncorrelated sequence
  ([NIST/SEMATECH e-Handbook, §1.3.3.1](@cite guthrie-2020-nist-ehandbook)).
  Outside it, a noise type is missing from the model or at the wrong
  level.

!!! note "A long record resolves small departures"
    Both intervals narrow as ``1/\\sqrt{n}``, so on a long record the test
    fails a model that departs a little from the noise and still predicts
    well. A model with flicker states is one, since its Markov processes
    only approximate a flicker spectrum; [What the states buy](@ref) has
    measured cases. Compare `rms` with `claimed` in [`holdover`](@ref)
    before discarding the model.

# Keywords

- `confidence`: probability that each statistic falls inside its interval
  under the model.

# Returns

A named tuple with

- `nis`: mean of ``e_k^2``, the normalized innovation squared;
- `nis_bounds`: its interval at `confidence`;
- `acf`: lag-1 autocorrelation of ``e_k``;
- `acf_bound`: the bound on its magnitude at `confidence`;
- `ok`: whether both are inside.

# Examples

```julia
using SigmaTau, Random
m = ClockModel(; tau0 = 1.0, R = 1e-22, q1 = 1e-24, q2 = 1e-30)
p, _ = simulate_clock(m, 4096; rng = Xoshiro(1))
innovation_check(kalman_filter(p, m))        # the model that made the record
innovation_check(kalman_filter(p, ClockModel(; tau0 = 1.0, R = 1e-22, q1 = 1e-25)))
```
"""
function innovation_check(kf::KalmanResult; confidence::Real=0.95)
    keep = findall(isfinite, kf.innovation)
    e = kf.innovation[keep] ./ sqrt.(kf.innovation_var[keep])
    n = length(e)
    n > 2 || throw(ArgumentError("innovation_check: needs more than two innovations"))
    χ² = Chisq(n)
    α = (1 - confidence) / 2
    bounds = (quantile(χ², α) / n, quantile(χ², 1 - α) / n)
    nis = sum(abs2, e) / n
    ē = e .- sum(e) / n
    acf = sum(ē[1:end-1] .* ē[2:end]) / sum(abs2, ē)
    bound = quantile(Normal(), 1 - α) / sqrt(n)
    return (nis=nis, nis_bounds=bounds, acf=acf, acf_bound=bound,
            ok=bounds[1] <= nis <= bounds[2] && abs(acf) <= bound)
end

"""
    simulate_clock(m::ClockModel, N; rng=nothing) → (PhaseData, Matrix{Float64})

A phase record of `N` samples drawn from clock model `m`, with the true
states behind it, for checking a filter against a known truth.

The clock starts at zero phase, frequency and drift, with its Markov states
drawn from their stationary distribution, and each step adds noise of
covariance [`process_noise`](@ref). A reading is the phase, plus the
Markov phase states, plus white noise of the model's reading variance.

# Arguments

- `m`: the clock model.
- `N`: number of samples.

# Keywords

- `rng`: a random number generator, such as `Xoshiro(1)`, for a
  reproducible record.

# Returns

The record, and the true states with one row per sample and the columns of
`x` in [`KalmanResult`](@ref).

# Examples

```julia
using SigmaTau, Random
m = ClockModel(; tau0 = 1.0, R = 1e-22, q1 = 1e-24, q2 = 1e-30)
p, X = simulate_clock(m, 4096; rng = Xoshiro(1))
kf = kalman_filter(p, m)
sqrt(sum(abs2, kf.x[:, 1] .- X[:, 1]) / 4096)    # RMS error of the phase estimate
```
"""
function simulate_clock(m::ClockModel, N::Int; rng=nothing)
    N >= 1 || throw(ArgumentError("simulate_clock: N must be ≥ 1, got $N"))
    np, n = _npoly(m), _nstates(m)
    Φ, Q = state_transition(m, m.tau0), process_noise(m, m.tau0)
    L = _psd_sqrt(Q)
    draw(k) = rng === nothing ? randn(k) : randn(rng, k)
    X = zeros(N, n)
    n > np && (X[1, np+1:n] = sqrt.(_markov_vars(m)) .* draw(n - np))
    for k in 2:N
        X[k, :] = Φ * X[k-1, :] .+ L * draw(n)
    end
    z = X * _reading_row(m) .+ sqrt(_reading_var(m)) .* draw(N)
    return PhaseData(z, m.tau0), X
end

# A square root of a covariance that may be singular (a noise left at zero).
function _psd_sqrt(Q::AbstractMatrix)
    F = eigen(Symmetric(Q))
    return F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
end

"""
    prediction_error(kf::KalmanResult, data::PhaseData, horizons; burnin=0)
        → Vector{Float64}

RMS error, in seconds, of the filter's predictions of the readings
`horizons` samples ahead.

From the estimate at each sample ``k`` past `burnin`, the state is carried
to ``k + h`` by [`state_transition`](@ref), and the reading it predicts is
compared with the reading there. The error therefore holds the noise of
that reading as well as the clock's. Pairs with a gap at either end are
skipped.

# Arguments

- `kf`: the filter's run over `data`; not a smoothed one, whose estimates
  already use the later readings.
- `data`: the record the filter ran on.
- `horizons`: prediction horizons in samples.

# Keywords

- `burnin`: samples at the start left out of the error.

# Examples

```julia
using SigmaTau, Random
m = ClockModel(; tau0 = 1.0, R = 1e-22, q1 = 1e-24, q2 = 1e-30)
p, _ = simulate_clock(m, 20_000; rng = Xoshiro(1))
prediction_error(kalman_filter(p, m), p, [1, 10, 100]; burnin = 1000)
```
"""
function prediction_error(kf::KalmanResult, data::PhaseData, horizons::AbstractVector{<:Integer}; burnin::Integer=0)
    z, X, n = data.x, kf.x, size(kf.x, 2)
    return map(horizons) do h
        w = state_transition(kf.model, h * kf.model.tau0)' * _reading_row(kf.model)
        s, c = 0.0, 0
        for k in burnin+1:length(z)-h
            isfinite(z[k+h]) && isfinite(z[k]) || continue
            s += abs2(z[k+h] - sum(X[k, j] * w[j] for j in 1:n)); c += 1
        end
        c > 0 || throw(ArgumentError("prediction_error: no readings $h samples past the burn-in"))
        sqrt(s / c)
    end
end

prediction_error(::KalmanResult, ::FrequencyData, ::AbstractVector{<:Integer}; kw...) = _needs_phase("prediction_error")

"""
    tune_q(data::PhaseData, m::ClockModel; burnin=length(data.x) ÷ 10,
           maxiter=300) → NamedTuple

Retune the white FM and random-walk FM levels of `m` by maximum likelihood.

The levels `q1` and `q2` are searched for the smallest
``-2\\ln L = \\sum_t [\\ln S_t + \\nu_t^2 / S_t]``, the sum over the
innovations ``\\nu_t`` of the filter and their variances ``S_t`` past the
burn-in ([Tryon and Jones (1983), Eq. (5.1)](@cite tryon-1983-parameters)),
with the other levels of `m` held fixed. Every trial is one pass of the
filter, and a missing reading adds no term. The search is Nelder–Mead on
the logarithms of the two levels, started from those of `m` with a step of
one natural logarithm, so it finds the optimum from a start whose ratio
`q1 / q2` is off by a factor of ten or more.

A noise fit gives the levels that match ``σ(τ)``; the likelihood gives the
levels under which the innovations have the variance the filter claims, and so
the error the filter claims one sample ahead is the error it makes. The two
differ when the record has noise that the model lacks, such as flicker FM.
The error the filter claims is then right at one sample and not beyond,
and in simulation it was wrong by a factor of up to 8 at longer holdover
times ([What the states buy](@ref)). A model with flicker states, the
default of [`clock_model`](@ref), carries that noise instead, and predicted
better there. See also [`innovation_check`](@ref).

# Arguments

- `data`: the phase record.
- `m`: the model to start from, with `q1 > 0` and `q2 > 0`.

# Keywords

- `burnin`: samples at the start left out of the sum.
- `maxiter`: iterations of the search at most.

# Returns

A named tuple `(model, cost, converged, evaluations)`: the tuned model;
``-2\\ln L`` at its levels; whether the search ended with the simplex and
the cost both within tolerance (otherwise `maxiter` was reached and
a restart from `model` may go on); and the number of filter runs.

# Examples

```julia
using SigmaTau, Random
h = Dict(0 => 2e-24, -1 => 1e-27, -2 => 1e-32)
p = noise_gen(PhaseData, 20_000, 1.0; h, rng = Xoshiro(1))
m = clock_model(powerlaw_fit(p; alphas = [2, 0, -2]))    # without flicker FM
t = tune_q(p, m)
t.model.q1 / m.q1, t.model.q2 / m.q2, t.converged
```
"""
function tune_q(data::PhaseData, m::ClockModel; burnin::Integer=length(data.x) ÷ 10, maxiter::Integer=300)
    m.q1 > 0 && m.q2 > 0 || throw(ArgumentError("tune_q: needs a model with q1 > 0 and q2 > 0 to start from"))
    length(data.x) > burnin || throw(ArgumentError("tune_q: the burn-in leaves no readings"))
    maxiter >= 1 || throw(ArgumentError("tune_q: maxiter must be at least 1, got $maxiter"))
    build(θ) = _with(m; q1=exp(θ[1]), q2=exp(θ[2]))
    function cost(θ)
        kf = kalman_filter(data, build(θ))
        s, n = 0.0, 0
        for k in burnin+1:length(data.x)
            ν, S = kf.innovation[k], kf.innovation_var[k]
            (isfinite(ν) && S > 0) || continue
            s += log(S) + ν^2 / S; n += 1
        end
        return n == 0 ? Inf : s
    end
    θ, c, converged, evaluations = _nelder_mead(cost, [log(m.q1), log(m.q2)]; maxiter)
    return (model=build(θ), cost=c, converged=converged, evaluations=evaluations)
end

tune_q(::FrequencyData, ::ClockModel; kw...) = _needs_phase("tune_q")

# Nelder–Mead with the standard coefficients (reflect 1, expand 2, contract
# ½, shrink ½), from `x0` and a step of `step` along each axis. It stops when
# the cost differs by less than `ftol` across the simplex and the simplex
# is smaller than `xtol`.
function _nelder_mead(f, x0::Vector{Float64}; step=1.0, maxiter=300, ftol=1e-3, xtol=1e-3)
    n = length(x0)
    P = [copy(x0) for _ in 1:n+1]
    for i in 1:n
        P[i+1][i] += step
    end
    F = f.(P)
    evals = n + 1
    converged = false
    for _ in 1:maxiter
        o = sortperm(F); P, F = P[o], F[o]
        if F[end] - F[1] < ftol && maximum(maximum(abs, P[i] - P[1]) for i in 2:n+1) < xtol
            converged = true
            break
        end
        c = sum(P[1:n]) / n
        xr = c + (c - P[end]); fr = f(xr); evals += 1
        if fr < F[1]
            xe = c + 2 * (c - P[end]); fe = f(xe); evals += 1
            P[end], F[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < F[n]
            P[end], F[end] = xr, fr
        else
            xc = fr < F[end] ? c + (xr - c) / 2 : c + (P[end] - c) / 2
            fc = f(xc); evals += 1
            if fc < min(fr, F[end])
                P[end], F[end] = xc, fc
            else
                for i in 2:n+1
                    P[i] = P[1] + (P[i] - P[1]) / 2; F[i] = f(P[i]); evals += 1
                end
            end
        end
    end
    i = argmin(F)
    return P[i], F[i], converged, evals
end
