# Single-clock Kalman filter on a phase record: filter, Rauch–Tung–Striebel
# smoother, holdover error and an innovation check, on a ClockModel.

"""
    KalmanResult

Output of [`kalman_filter`](@ref) or [`kalman_smooth`](@ref): the state
estimate at every sample, its covariance, and the innovations.

$(TYPEDFIELDS)
"""
struct KalmanResult
    "Clock model the filter ran on."
    model::ClockModel
    "Sample times in seconds from the first sample."
    t::Vector{Float64}
    "State estimates, one row per sample: phase (s), frequency, and drift (1/s) with a drift state."
    x::Matrix{Float64}
    "State covariances, `P[:, :, k]` at sample k."
    P::Array{Float64,3}
    "Innovations z − H x̂⁻ (s); `NaN` at gaps."
    innovation::Vector{Float64}
    "Innovation variances H P⁻ Hᵀ + R (s²); `NaN` at gaps."
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

"""
    kalman_filter(data::PhaseData, m::ClockModel; x0, P0) → KalmanResult

Run the Kalman filter for clock model `m` over a phase record. The state is
phase and frequency, plus drift when `m.drift`; each sample is a phase
reading with white PM variance `m.R`. `NaN` samples are gaps: the filter
predicts across them without an update.

The start `x0`, `P0` defaults to a least-squares polynomial through the
first samples, with its covariance inflated tenfold. The covariance update
uses the Joseph form, which keeps P symmetric and positive under round-off.

Pass the result to [`kalman_smooth`](@ref), [`holdover_error`](@ref) or
[`innovation_check`](@ref).

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 4096, 1.0; h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
m = clock_model(powerlaw_fit(p; alphas = [2, 0, -2]))
kf = kalman_filter(p, m)
kf.x[end, 2]               # frequency estimate at the last sample
holdover_error(kf, [60.0, 3600.0])
```
"""
function kalman_filter(data::PhaseData, m::ClockModel;
                       x0::Union{AbstractVector,Nothing}=nothing,
                       P0::Union{AbstractMatrix,Nothing}=nothing)
    isapprox(data.tau0, m.tau0; rtol=1e-9) || throw(ArgumentError(
        "kalman_filter: record τ₀ = $(data.tau0) s but the model's is $(m.tau0) s"))
    z = Float64.(data.x)
    N, n = length(z), m.drift ? 3 : 2
    xs, Ps0 = x0 === nothing || P0 === nothing ? _kalman_start(z, m.tau0, n) : (x0, P0)
    x̂ = Float64.(x0 === nothing ? xs : x0)
    P = Float64.(P0 === nothing ? Ps0 : P0)
    size(P) == (n, n) && length(x̂) == n || throw(ArgumentError("kalman_filter: x0 and P0 must have $n states"))
    Φ, Q = state_transition(m, m.tau0), process_noise(m, m.tau0)
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
            S[k] = P[1, 1] + m.R
            ν[k] = z[k] - x̂[1]
            K = P[:, 1] ./ S[k]
            x̂ = x̂ .+ K .* ν[k]
            IKH = Matrix{Float64}(I, n, n); IKH[:, 1] .-= K
            P = IKH * P * IKH' .+ m.R .* (K * K')
            P = (P + P') ./ 2
        end
        X[k, :] = x̂; Ps[:, :, k] = P
    end
    t = (0:N-1) .* m.tau0
    return KalmanResult(m, collect(t), X, Ps, ν, S, false, Xp, Pps)
end

"""
    kalman_smooth(kf::KalmanResult) → KalmanResult

Rauch–Tung–Striebel smoother over a filtered run: each state estimate then
uses the whole record, before and after it. The innovations are the
filter's.
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

Standard deviation (s) of the time error after free-running for each τ in
`tau` (seconds) from the last sample, with the filter's last estimate as the
starting point: `√[Φ(τ) P Φ(τ)ᵀ + Q(τ)]₁₁`.
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

Whether the model fits the record, from the normalized innovations
eₖ = νₖ/√Sₖ, which are white with unit variance when it does:

- `nis`: mean of eₖ², with the χ² interval `nis_bounds` it falls in at
  `confidence` when the model is right. Above the interval, the model's
  noise is too low; below, too high.
- `acf`: lag-1 autocorrelation of eₖ, with `acf_bound`, the ±bound it
  stays within for white innovations. Outside it, a noise type is missing
  or mislevelled.
- `ok`: both inside.
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
    simulate_clock(m::ClockModel, N; rng=nothing) → (PhaseData, Matrix)

A phase record of `N` samples from clock model `m`, and its true states
(one row per sample), for checking filters against a known truth. The
clock starts at zero phase, frequency and drift.
"""
function simulate_clock(m::ClockModel, N::Int; rng=nothing)
    n = m.drift ? 3 : 2
    Φ, Q = state_transition(m, m.tau0), process_noise(m, m.tau0)
    L = _psd_sqrt(Q)
    draw(k) = rng === nothing ? randn(k) : randn(rng, k)
    X = zeros(N, n)
    for k in 2:N
        X[k, :] = Φ * X[k-1, :] .+ L * draw(n)
    end
    z = X[:, 1] .+ sqrt(m.R) .* draw(N)
    return PhaseData(z, m.tau0), X
end

# A square root of a covariance that may be singular (a noise left at zero).
function _psd_sqrt(Q::AbstractMatrix)
    F = eigen(Symmetric(Q))
    return F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
end

"""
    prediction_error(kf::KalmanResult, data::PhaseData, horizons; burnin=0) → Vector{Float64}

RMS error (s) of the filter's phase predictions `horizons` samples ahead:
from the estimate at each sample k past `burnin`, extrapolate the phase to
k + h and compare it with the reading there. Gaps are skipped.
"""
function prediction_error(kf::KalmanResult, data::PhaseData, horizons::AbstractVector{<:Integer}; burnin::Integer=0)
    z, X, n = data.x, kf.x, size(kf.x, 2)
    return map(horizons) do h
        δ = h * kf.model.tau0
        w = [δ^(j - 1) / factorial(j - 1) for j in 1:n]
        s, c = 0.0, 0
        for k in burnin+1:length(z)-h
            isfinite(z[k+h]) && isfinite(z[k]) || continue
            s += abs2(z[k+h] - sum(X[k, j] * w[j] for j in 1:n)); c += 1
        end
        c > 0 || throw(ArgumentError("prediction_error: no readings $h samples past the burn-in"))
        sqrt(s / c)
    end
end

"""
    tune_q(data::PhaseData, m::ClockModel; horizons=[10, 100, 1000], decades=2, points=9, burnin)
        → (model, q1, q2, cost)

Retune the white FM and random-walk FM levels of `m` for prediction: run
the filter over a log grid of `points` values spanning ± `decades` around
`m.q1` and `m.q2`, and keep the pair whose phase predictions
`horizons` samples ahead have the smallest summed log RMS error (see
[`prediction_error`](@ref)). `R` and `q3` stay as they are. The burn-in
defaults to a tenth of the record. Returns the tuned model, the grid axes
and the cost on the grid (rows q1, columns q2).

A noise fit gives levels that match σ(τ); tuning gives the ones that
predict best, which differ when the record has noise the model lacks, such
as flicker FM.
"""
function tune_q(data::PhaseData, m::ClockModel; horizons::AbstractVector{<:Integer}=[10, 100, 1000],
                decades::Real=2, points::Integer=9, burnin::Integer=length(data.x) ÷ 10)
    m.q1 > 0 && m.q2 > 0 || throw(ArgumentError("tune_q: needs a model with q1 > 0 and q2 > 0 to centre the grid"))
    hs = filter(h -> h < length(data.x) - burnin, horizons)
    isempty(hs) && throw(ArgumentError("tune_q: every horizon is longer than the record past the burn-in"))
    grid(q) = q .* exp10.(range(-decades, decades; length=points))
    q1s, q2s = grid(m.q1), grid(m.q2)
    cost = [begin
                mi = ClockModel(m.tau0, m.R, a, b, m.q3, m.drift)
                sum(log, prediction_error(kalman_filter(data, mi), data, hs; burnin))
            end for a in q1s, b in q2s]
    i = argmin(cost)
    return (model=ClockModel(m.tau0, m.R, q1s[i[1]], q2s[i[2]], m.q3, m.drift), q1=q1s, q2=q2s, cost=cost)
end
