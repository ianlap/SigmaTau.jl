# Clock steering: frequency steers computed from the Kalman estimate of a
# clock's phase and frequency against its reference, with the gains set by
# where they put the poles of the closed loop (Matsakis 2019; Matsakis 2020
# for the integral gain). A measured record is replayed through
# the loop; a clock model in its place gives the steady state it settles to.

"""
    SteeringGains(phase, frequency, integral, interval)

Gains of a steering loop. Every `interval` seconds the frequency of the
clock is changed by

```math
u = -\\left(g_P\\,\\hat p + g_D\\,\\hat y + g_I \\sum \\hat p\\right),
```

with ``\\hat p`` and ``\\hat y`` the Kalman estimates of its phase and
frequency against the reference and ``\\sum \\hat p`` the sum of the phase
estimates at the steers so far, this one included. The change stays in the
clock until the next steer adds to it. ``g_P`` and ``g_D`` are the gains
``g_1`` and ``g_2`` of
[Matsakis (2019), Eqs. (A.17)–(A.18)](@cite matsakis-2019-steering); the
integral gain ``g_I`` is that of
[Matsakis (2020), Eq. (3.1)](@cite matsakis-2020-pid). Build the gains with
[`steering_gains`](@ref) or [`lqg_gains`](@ref).

# Fields

$(TYPEDFIELDS)
"""
struct SteeringGains
    "Phase gain (1/s): the frequency steer per second of phase error."
    phase::Float64
    "Frequency gain: the share of the frequency error one steer takes out."
    frequency::Float64
    "Integral gain (1/s), on the running sum of the phase error; `0` for a loop of two gains."
    integral::Float64
    "Seconds between steers."
    interval::Float64
    function SteeringGains(phase, frequency, integral, interval)
        interval > 0 || throw(ArgumentError("SteeringGains: interval must be positive, got $interval"))
        return new(phase, frequency, integral, interval)
    end
end

Base.show(io::IO, g::SteeringGains) =
    print(io, "SteeringGains(phase=", _show4g(g.phase), " /s, frequency=", _show4g(g.frequency),
          g.integral == 0 ? "" : ", integral=" * _show4g(g.integral) * " /s", ", every ", _show4g(g.interval), " s)")

"""
    steering_gains(T, interval; integral=false) → SteeringGains
    steering_gains(poles::AbstractVector, interval) → SteeringGains

Gains of a steering loop by pole placement.

With a time constant `T` the loop is critically damped: its poles sit
together at ``r = e^{-τ/T}``, with τ the interval between steers, and a
disturbance dies away as ``e^{-t/T}`` without ringing. The gains are those
of [Matsakis (2019), Eqs. (3.10)–(3.11)](@cite matsakis-2019-steering),

```math
g_P = (1 - r)^2/τ, \\qquad g_D = 1 - r^2 ,
```

and with `integral = true`, all three poles at ``r``, those of
[Matsakis (2020), Eqs. (3.15)–(3.17)](@cite matsakis-2020-pid),

```math
g_I = (1 - r)^3/τ, \\qquad g_P = (1 - 3r^2 + 2r^3)/τ, \\qquad g_D = 1 - r^3 .
```

`poles` places them anywhere instead: two poles for phase and frequency
gains, three with an integral gain. With ``w_i = 1 - z_i`` for each pole
``z_i``, the gains are ``τ g_P = w_1 w_2`` and ``g_D = w_1 + w_2 - w_1 w_2``
for two poles, and for three

```math
τ g_I = w_1 w_2 w_3, \\qquad τ g_P = e_2 - 2 w_1 w_2 w_3, \\qquad
g_D = e_1 - e_2 + w_1 w_2 w_3 ,
```

with ``e_1 = \\sum w_i`` and ``e_2 = \\sum_{i<j} w_i w_j``. These follow from
the characteristic polynomial given under [`steering_poles`](@ref).

# Arguments

- `T`: time constant of the loop in seconds.
- `interval`: seconds between steers.
- `poles`: two or three poles in the z-plane; complex poles in conjugate pairs.

# Keywords

- `integral`: whether the loop has an integral gain.

# Examples

```jldoctest
julia> using SigmaTau

julia> g = steering_gains(3600.0, 60.0)    # settle in an hour, steer each minute
SteeringGains(phase=4.553e-06 /s, frequency=0.03278, every 60 s)

julia> round.(steering_poles(g).time_constants)
2-element Vector{Float64}:
 3600.0
 3600.0
```
"""
function steering_gains(T::Real, interval::Real; integral::Bool=false)
    T > 0 && interval > 0 || throw(ArgumentError("steering_gains: T and interval must be positive, got $T and $interval"))
    return _gains_from(fill(complex(-expm1(-interval / T)), integral ? 3 : 2), interval)
end

function steering_gains(poles::AbstractVector{<:Number}, interval::Real)
    length(poles) in (2, 3) || throw(ArgumentError(
        "steering_gains: two poles for phase and frequency gains, three with an integral gain; got $(length(poles))"))
    return _gains_from(ComplexF64[1 - z for z in poles], interval)
end

# Gains from w = 1 − z; a third w of zero is the integrator left out of the loop.
function _gains_from(w::Vector{ComplexF64}, interval::Real)
    w₃ = length(w) == 3 ? w[3] : zero(ComplexF64)
    e₁, e₂, e₃ = w[1] + w[2] + w₃, w[1] * w[2] + (w[1] + w[2]) * w₃, w[1] * w[2] * w₃
    g = (e₂ - 2e₃, e₁ - e₂ + e₃, e₃)
    all(x -> abs(imag(x)) <= 1e-9 * max(abs(x), abs(e₁)), g) || throw(ArgumentError(
        "steering_gains: complex poles must come in conjugate pairs"))
    return SteeringGains(real(g[1]) / interval, real(g[2]), real(g[3]) / interval, interval)
end

"""
    lqg_gains(interval; phase, frequency=0, steer) → SteeringGains

Gains of the steering loop that minimizes a quadratic cost: the linear
quadratic Gaussian (LQG) design of
[Koppang and Leland (1999), §III](@cite koppang-1999-lqg-masers).

The cost is the sum over the steers of

```math
W_p\\,p^2 + W_f\\,f^2 + W_u\\,u^2 ,
```

with ``p`` and ``f`` the phase and frequency of the clock against its
reference and ``u`` the steer. The gains are
``G = (B^T K B + W_u)^{-1} B^T K Φ``, with ``K`` the solution of the Riccati
equation

```math
K = Φ^T K Φ + W_Q - Φ^T K B\\,(B^T K B + W_u)^{-1} B^T K Φ, \\qquad
W_Q = \\begin{pmatrix} W_p & 0 \\\\ 0 & W_f \\end{pmatrix},
```

and ``Φ`` and ``B`` those of the loop at the interval τ between steers.
They do not depend on the noise of the clock or of the readings, and only
the ratios ``W_p τ^2/W_u`` and ``W_f/W_u`` matter
([Matsakis (2019), §1](@cite matsakis-2019-steering)). A small steer weight
gives ``g_P = 1/τ`` and ``g_D = 1``, which takes out the whole estimated
error in one steer; a large one gives a slow loop. The design has no
integral gain: the cost is least without one
([Matsakis (2020), §III](@cite matsakis-2020-pid)).

# Arguments

- `interval`: seconds between steers.

# Keywords

- `phase`: weight ``W_p`` on the square of the phase (s), positive.
- `frequency`: weight ``W_f`` on the square of the frequency.
- `steer`: weight ``W_u`` on the square of the steer, positive.

# Examples

The first design of [Koppang and Leland (1999), §IV](@cite koppang-1999-lqg-masers),
for a steer each hour, which prints the gains as `6.503e-5` and `0.5771`:

```jldoctest
julia> using SigmaTau

julia> lqg_gains(3600.0; phase = 1e-3, frequency = 1e-3, steer = 1e5)
SteeringGains(phase=6.503e-05 /s, frequency=0.5771, every 3600 s)
```
"""
function lqg_gains(interval::Real; phase::Real, frequency::Real=0, steer::Real)
    interval > 0 && phase > 0 && frequency >= 0 && steer > 0 || throw(ArgumentError(
        "lqg_gains: interval, phase and steer must be positive and frequency not negative"))
    # In units of the interval, with the steer weight as the unit of cost:
    # state (p, τ f), steer τ u, Φ = [1 1; 0 1], B = (1, 1).
    Φ, b = [1.0 1.0; 0.0 1.0], [1.0, 1.0]
    WQ = [phase*interval^2/steer 0.0; 0.0 frequency/steer]
    K = _riccati(Φ, b * b', WQ)
    g = (Φ' * K * b) ./ (dot(b, K, b) + 1)
    return SteeringGains(g[1] / interval, g[2], 0.0, interval)
end

"""
    steering_poles(g::SteeringGains) → NamedTuple

Where the gains `g` put the poles of the closed loop, and what follows from
them. The poles are the roots of

```math
(z - 1)^3 + τ g_I\\,z^2 + τ g_P\\,z (z - 1) + g_D\\,(z - 1)^2 = 0 ,
```

with τ the interval between steers
([Matsakis (2020), Eqs. (3.10) and (4.15)](@cite matsakis-2020-pid)).
Without an integral gain a factor ``z - 1`` drops out and leaves the
quadratic of
[Matsakis (2019), Eqs. (3.8)–(3.9)](@cite matsakis-2019-steering).

# Returns

A named tuple with

- `poles`: the roots, complex;
- `time_constants`: for each pole the ``T`` (s) of its decay ``e^{-t/T}``,
  from ``|z| = e^{-τ/T}``, negative for a pole that grows;
- `frequencies`: the frequency (Hz) each pole rings at, zero on the positive
  real axis and ``1/2τ`` on the negative;
- `stable`: whether every pole is inside the unit circle.

For two gains the loop is stable when ``g_P > 0``, ``g_D > 0`` and
``τ g_P + 2 g_D < 4``
([Matsakis (2019), Eq. (3.17) and the text after Eq. (3.18)](@cite matsakis-2019-steering)).

!!! note "Repeated poles"
    A pole placed twice or three times is found to about the square or cube
    root of machine precision, as a fraction of its distance from 1.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = steering_poles(SteeringGains(0.02, 0.3, 0, 2.0));

julia> p.stable, round.(p.time_constants; digits = 3)
(true, [11.215, 11.215])
```
"""
function steering_poles(g::SteeringGains)
    τ = g.interval
    P, D, J = τ * g.phase, g.frequency, τ * g.integral
    # In w = 1 − z the polynomial is w² − (P + D) w + P, or with an integral
    # gain J: w³ − (J + P + D) w² + (2J + P) w − J. Its companion matrix:
    n = J == 0 ? 2 : 3
    C = zeros(ComplexF64, n, n)
    for i in 2:n
        C[i, i-1] = 1
    end
    if n == 2
        C[1, 2], C[2, 2] = -P, P + D
    else
        C[1, 3], C[2, 3], C[3, 3] = J, -(2J + P), J + P + D
    end
    w = Vector{ComplexF64}(eigvals(C; sortby=nothing))
    logz = [log1p(-2real(x) + abs2(x)) / 2 for x in w]           # log |1 − w|
    return (poles=1 .- w, time_constants=-τ ./ logz,
            frequencies=[abs(angle(1 - x)) / (2π * τ) for x in w], stable=all(<(0), logz))
end

# What the gains act on, as rows on the filter's state at a steer every τ
# seconds: the frequency ŷ, and the steer made ahead of any feedback. ŷ is
# the frequency state plus the Markov frequency states. A drift state d is
# steered out ahead, u = … − τ d̂, and ŷ taken at mid-interval less that
# steer, f̂ − τ d̂/2: phase and ŷ then obey the two-state loop exactly, and a
# steady drift leaves no standing offset.
function _steering_rows(m::ClockModel, τ::Float64)
    n, np = _nstates(m), _npoly(m)
    y, ahead = zeros(n), zeros(n)
    y[2] = 1.0
    y[np+1:np+length(_fm_taus(m))] .= 1.0
    m.drift && (y[3] = -τ / 2; ahead[3] = τ)
    return y, ahead
end

"""
    SteeringResult

Output of [`steer`](@ref) on a record: the record as the loop would have
left it, and what the loop did.

# Fields

$(TYPEDFIELDS)
"""
struct SteeringResult
    "Gains of the loop."
    gains::SteeringGains
    "The steered record: each reading plus the phase the steers had added by then."
    steered::PhaseData{Float64}
    "Frequency steer made at each sample; zero between steers."
    steers::Vector{Float64}
    "Sum of the steers so far: the frequency correction in force at each sample."
    correction::Vector{Float64}
    "The filter's estimate of the steered clock's phase (s) at each sample."
    phase::Vector{Float64}
    "The filter's estimate of the steered clock's frequency at each sample."
    frequency::Vector{Float64}
end

Base.show(io::IO, r::SteeringResult) =
    print(io, "SteeringResult(", length(r.steers), " samples, ", r.gains, ")")

"""
    steer(data::PhaseData, model::ClockModel, g::SteeringGains) → SteeringResult

Replay a record through a steering loop: what the clock would have done had
it been steered to its reference with the gains `g`.

`data` is the phase of the free-running clock against the reference. The
Kalman filter for `model` follows the clock ([`kalman_filter`](@ref)), and
every `g.interval` the frequency of the clock is changed by the ``u`` of
[`SteeringGains`](@ref), computed from the estimates just after the reading
and applied at once, as in
[Koppang and Leland (1999), Eqs. (1)–(3)](@cite koppang-1999-lqg-masers) and
[Matsakis (2019), Eqs. (A.11)–(A.17)](@cite matsakis-2019-steering). The
steered record is each reading plus the phase the steers have added by
then, so any deviation of it shows the stability of the steered clock.

The filter runs at every sample and the steers come every `g.interval`, a
whole multiple of the sample interval, with gains computed for that
interval ([Gödel et al. (2019), §3.1](@cite goedel-2019-steering-comparison)).
The filter knows each steer exactly, so its estimate of the steered clock
is its estimate of the free-running one plus the correction, and one pass
of the filter serves any gains.

!!! note "Beyond phase and frequency"
    The published loops act on a phase and a frequency state. Here the
    frequency ``\\hat y`` is the whole frequency estimate of the model: the
    frequency state plus the Markov states of a model with flicker FM. A
    drift state ``\\hat d`` is steered out ahead of the feedback,
    ``u = … - τ \\hat d`` with ``\\hat y`` taken as ``\\hat f - τ \\hat d/2``,
    which leaves phase and ``\\hat y`` obeying the two-state loop, so a steady
    drift leaves no standing phase offset. Both extensions are this
    package's.

!!! warning "A drift the model does not hold"
    Without a drift state the filter lags a drifting clock, and the steered
    phase stands off from the reference. An integral gain brings the mean of
    the estimated phase to zero, not that lag.

# Arguments

- `data`: phase of the free-running clock against its reference; `NaN`
  samples are gaps, which the filter predicts across.
- `model`: clock model of the pair, with `model.tau0` the sample interval.
- `g`: gains, with `g.interval` a whole multiple of `model.tau0`.

# Examples

```julia
using SigmaTau, Random
m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26)
p, _ = simulate_clock(m, 65536; rng = Xoshiro(1))
r = steer(p, m, steering_gains(300.0, 1.0))
adev(r.steered)             # the steered clock against its reference
```
"""
function steer(data::PhaseData, m::ClockModel, g::SteeringGains)
    τ₀ = m.tau0
    n = round(Int, g.interval / τ₀)
    n >= 1 && isapprox(n * τ₀, g.interval; rtol=1e-9) || throw(ArgumentError(
        "steer: the gains' interval of $(g.interval) s is not a whole multiple of tau0 = $τ₀ s"))
    kf = kalman_filter(data, m)
    y, ahead = _steering_rows(m, g.interval)
    N = length(kf.t)
    c, s, u, p̂, ŷ = zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)
    Σp = 0.0
    for k in 1:N
        if k > 1
            s[k] = s[k-1] + u[k-1]
            c[k] = c[k-1] + τ₀ * s[k]
        end
        x̂ = view(kf.x, k, :)
        p̂[k], ŷ[k] = x̂[1] + c[k], dot(y, x̂) + s[k]
        k % n == 0 || continue
        Σp += p̂[k]
        u[k] = -(g.phase * p̂[k] + g.frequency * ŷ[k] + g.integral * Σp + dot(ahead, x̂))
    end
    steered = PhaseData(Float64.(data.x) .+ c, τ₀; source=data.source, start=data.start)
    return SteeringResult(g, steered, u, s, p̂, ŷ)
end

"""
    SteeringPrediction

Output of [`steer`](@ref) on a clock model: the steady state of the loop.

# Fields

$(TYPEDFIELDS)
"""
struct SteeringPrediction
    "Gains of the loop."
    gains::SteeringGains
    "RMS time error (s) of the steered clock against its reference."
    phase::Float64
    "RMS frequency error of the steered clock; its white FM, which has no value at an instant, is left out."
    frequency::Float64
    "RMS frequency steer."
    steers::Float64
    "Averaging times (s)."
    tau::Vector{Float64}
    "Allan deviation of the steered clock against its reference at `tau`."
    adev::Vector{Float64}
end

Base.show(io::IO, p::SteeringPrediction) =
    print(io, "SteeringPrediction(phase ", _show4g(p.phase), " s, frequency ", _show4g(p.frequency),
          ", steers ", _show4g(p.steers), " rms; ", p.gains, ")")

# Σ = M Σ Mᵀ + W for a stable M, as the sum of Mᵏ W Mᵏᵀ: each pass doubles
# the number of terms.
function _lyapunov(M::Matrix{Float64}, W::Matrix{Float64})
    Σ, A = copy(W), copy(M)
    for _ in 1:128
        Δ = A * Σ * A'
        Σ .+= Δ
        all(abs(Δ[i]) <= eps() * abs(Σ[i]) for i in eachindex(Δ)) && break
        A = A * A
    end
    return (Σ + Σ') ./ 2
end

# X = Aᵀ X (I + G X)⁻¹ A + H, the discrete Riccati equation, by doubling:
# each pass stands for twice as many steps of the recursion as the last,
# which reaches the limit however slowly the recursion itself would.
function _riccati(A::Matrix{Float64}, G::Matrix{Float64}, H::Matrix{Float64})
    for _ in 1:128
        W = inv(I + G * H)
        H′ = H + A' * H * W * A
        G = G + A * W * G * A'
        A = A * W * A
        done = maximum(abs, H′ - H) <= 4eps() * maximum(abs, H′)
        H = (H′ + H′') ./ 2
        done && break
    end
    return H
end

# The Kalman filter of a model without drift once it has settled: the
# covariance before a reading, the gain, and the variance of the
# innovations. It works in units of τ₀ and of the phase σ one step adds, and
# a reading without noise is given a millionth of the clock's.
function _steady_filter(m::ClockModel)
    τ = m.tau0
    Φ, Q, h, R = state_transition(m, τ), process_noise(m, τ), _reading_row(m), _reading_var(m)
    n = length(h)
    σ = sqrt(R + Q[1, 1])
    σ > 0 || throw(ArgumentError("steer: the model has no noise"))
    d = fill(σ / τ, n)
    d[1] = σ
    d[n-length(_pm_taus(m))+1:n] .= σ
    P = _riccati(Matrix(((Φ .* d') ./ d)'), (h .* d) * (h .* d)' ./ max(R, 1e-12 * Q[1, 1]), Q ./ (d * d')) .* (d * d')
    S = dot(h, P, h) + R
    K = P * h ./ S
    m.q2 == 0 && (K[2] = 0.0)                    # a frequency that never walks is known by now
    return (Φ=Φ, Q=Q, h=h, R=R, P=P, K=K, S=S)
end

"""
    steer(model::ClockModel, g::SteeringGains; tau,
          readings=false) → SteeringPrediction

The steady state a steering loop settles to on a clock with the noise of
`model`: the RMS of the steered clock's time and frequency error and of the
steers, and its Allan deviation against the reference.

Steering does not change what the filter knows, so its gain ``K`` and its
covariances are those of the free-running clock, and its estimate obeys
``\\hat x_{k+1} = (Φ - B G)\\,\\hat x_k + K ν_{k+1}`` with white innovations
``ν``. The variances then follow from Lyapunov equations, as in
[Matsakis (2019), §2 and Appendix B](@cite matsakis-2019-steering) and, with
an integral gain, [Matsakis (2020), §II](@cite matsakis-2020-pid). The
figures are for the clock itself, without the noise of the readings that a
measured record carries; `readings = true` gives the time error and the
Allan deviation of the readings, to set beside a replayed record.

The loop steers at every sample here, `g.interval == model.tau0`; for a
longer interval replay a simulated record. A drift state is taken as
steered out. Random-run FM (`q3 > 0`) has no steady state.

# Arguments

- `model`: clock model of the pair of clock and reference.
- `g`: gains of a stable loop.

# Keywords

- `tau`: averaging times (s) of the Allan deviation, whole multiples of
  `model.tau0`; by default the octaves up to 32 times the longest time
  constant of the loop.
- `readings`: `phase` and `adev` are those of the readings of the steered
  clock, with their noise (white PM, and flicker PM where the model holds
  it), as in the `steered` record of `steer(data, model, g)`.

# Examples

```julia
using SigmaTau
m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26)
s = steer(m, steering_gains(300.0, 1.0))
s.phase                # RMS time error against the reference, s
s.adev                 # at s.tau
```
"""
function steer(m::ClockModel, g::SteeringGains; tau::Union{AbstractVector{<:Real},Nothing}=nothing, readings::Bool=false)
    τ = m.tau0
    isapprox(g.interval, τ; rtol=1e-9) || throw(ArgumentError(
        "steer: the steady state is worked out for a steer at every sample; the gains' interval is $(g.interval) s and tau0 is $τ s"))
    m.q3 == 0 || throw(ArgumentError("steer: random-run FM (q3 > 0) has no steady state"))
    poles = steering_poles(g)
    poles.stable || throw(ArgumentError("steer: these gains are unstable; there is no steady state"))
    m.drift && (m = ClockModel(τ, m.R, m.q1, m.q2, 0.0, false, m.flicker_fm, m.fm_band, m.flicker_pm, m.pm_band))
    (; Φ, Q, h, R, P, K) = _steady_filter(m)
    n = length(h)
    y, _ = _steering_rows(m, τ)
    gx = g.frequency .* y
    gx[1] += g.phase
    # ζ = (Σp̂, x̂, e), e the error of the estimate; without an integral gain
    # Σp̂ is left out. The noise is the clock's w and the reading's v:
    #   x̂⁺ = (Φ − b gxᵀ) x̂ − b g_I Σp̂ + K ν⁺,   ν⁺ = hᵀ(Φ e + w) + v
    #   e⁺ = (I − K hᵀ)(Φ e + w) − K v,          Σp̂⁺ = Σp̂ + p̂⁺
    o = g.integral == 0 ? 0 : 1
    xs, es = o+1:o+n, o+n+1:o+2n
    M, B = zeros(o + 2n, o + 2n), zeros(o + 2n, n + 1)
    b = Φ[:, 2]                                   # a steer moves the frequency, and the phase with it over the step
    M[xs, xs] = Φ .- b * gx'
    M[xs, es] = K * (h' * Φ)
    M[es, es] = (I - K * h') * Φ
    B[xs, 1:n], B[xs, n+1] = K * h', K
    B[es, 1:n], B[es, n+1] = I - K * h', -K
    if o == 1
        M[xs, 1] = -g.integral .* b
        M[1, :] = M[o+1, :]
        M[1, 1] += 1
        B[1, :] = B[o+1, :]
    end
    N = zeros(n + 1, n + 1)
    N[1:n, 1:n], N[n+1, n+1] = Q, R
    Σ = _lyapunov(M, B * N * B')
    cp, cy, cu = zeros(o + 2n), zeros(o + 2n), zeros(o + 2n)
    # The phase of the clock, or its reading z = hᵀx + v. The noise v of a
    # reading is not in the clock then, ⟨xₖ vₖ⟩ = 0, and the steers have put
    # it there j steps on: ⟨ζₖ₊ⱼ vₖ⟩ = Mʲ bv.
    if readings
        cp[xs], cp[es] = h, h
    else
        cp[o+1] = cp[o+n+1] = 1.0
    end
    Rz = readings ? R : 0.0
    bv = B[:, n+1] .* Rz
    cy[xs], cy[es] = y, y
    cu[xs] = -gx
    o == 1 && (cu[1] = -g.integral)
    ms = tau === nothing ? _grid(Octave, ceil(Int, 32 * maximum(poles.time_constants) / τ)) :
         [(k = round(Int, t / τ); k >= 1 && isapprox(k * τ, t; rtol=1e-9) ? k :
           throw(ArgumentError("steer: tau = $t is not a whole multiple of tau0 = $τ"))) for t in tau]
    # σ_y²(jτ) = ⟨(p₂ⱼ − 2pⱼ + p₀)²⟩ / 2(jτ)², from the covariance of the phase j and 2j steps apart.
    σy = map(ms) do j
        r₁ = (M^j)' * cp
        r₂ = (M^j)' * r₁
        sqrt(max(3dot(cp, Σ, cp) - 4dot(r₁, Σ, cp) + dot(r₂, Σ, cp) + 3Rz - 4dot(r₁, bv) + dot(r₂, bv), 0.0)) / (j * τ)
    end
    return SteeringPrediction(g, sqrt(dot(cp, Σ, cp) + Rz), sqrt(dot(cy, Σ, cy)), sqrt(dot(cu, Σ, cu)), ms .* τ, σy)
end
