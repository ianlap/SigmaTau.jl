# Time scales from a clock ensemble. EWFA, AT1 and AT2 average the clocks'
# predicted phases, with equal or learned inverse-variance weights; the Kalman
# scales run one filter over all the clocks, or one per clock difference.

"""
    TimescaleResult

Output of [`timescale`](@ref): each clock against the scale at every
sample, with its frequency, weight and prediction error.

[`timescale_phase`](@ref) gives the scale itself, against the reference of
the comparisons. Before a clock's first reading its `x`, `y` and `xhat`
are `NaN`.

# Fields

$(TYPEDFIELDS)
"""
struct TimescaleResult
    "Algorithm that formed the scale (`:ewfa`, `:at1`, `:at2`, `:greenhall`, `:jones_tryon`, `:stein`)."
    algorithm::Symbol
    "The clocks the scale was formed from."
    comparisons::ClockComparisons
    "Sample times in seconds from the first sample."
    t::Vector{Float64}
    "Each clock minus the scale (s), one column per clock."
    x::Matrix{Float64}
    "Each clock's frequency against the scale, as the algorithm estimated it."
    y::Matrix{Float64}
    "Weight of each clock at each sample: a row sums to 1, or is zero where no clock of the scale is read; `NaN` for `:jones_tryon`, which has no weights."
    weights::Matrix{Float64}
    "Each clock's prediction error (s); `NaN` where it had no reading."
    error::Matrix{Float64}
    "Each clock's phase against the scale as the algorithm estimated it (s): the filter's phase state for `:greenhall`, `:jones_tryon` and `:stein`, where `x` is read at the pivot; `x` itself for the others."
    xhat::Matrix{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", r::TimescaleResult)
    n, k = size(r.x)
    print(io, uppercase(String(r.algorithm)), " timescale from ", k, " clocks, ", n, " samples")
    @printf(io, "\n  final weights %s", join((@sprintf("%.3f", w) for w in r.weights[end, :]), ", "))
end

"""
    timescale_phase(r::TimescaleResult) → PhaseData

The reference of the comparisons against the scale, as a phase record for
the usual stability analysis.

At each sample it is ``x_i - z_i`` for the first clock ``i`` read then:
the clock against the scale, `r.x`, less the reading of the clock against
the reference. The clocks read at one sample give the same value, and a
sample at which none was read is `NaN`. Its deviations are those of the
scale as far as the reference is the more stable of the two. Against the
ideal reference of [`simulate_ensemble`](@ref) with `reference = 0` it is
the error of the scale, with the sign reversed.

# Examples

Two clocks with equal weights: the scale runs midway between them.

```jldoctest
julia> using SigmaTau

julia> cc = ClockComparisons([0.0 2.0; 1.0 5.0; 2.0 8.0; 3.0 11.0] .* 1e-9, 1.0);

julia> u = timescale_phase(timescale(cc; algorithm = :ewfa));

julia> u.x ≈ -[1.0, 3.0, 5.0, 7.0] .* 1e-9
true
```
"""
function timescale_phase(r::TimescaleResult)
    c = r.comparisons.x
    u = Float64[(i = findfirst(isfinite, view(c, n, :))) === nothing ? NaN : r.x[n, i] - c[n, i]
                for n in axes(c, 1)]
    return PhaseData(u, r.comparisons.tau0; source=uppercase(String(r.algorithm)))
end

"""
    implicit_mean(r::TimescaleResult, truth) → PhaseData

The implicit ensemble mean of a `:greenhall` scale against true time, for
a simulated ensemble whose true phases are known.

```math
\\bar{x}_0 = \\sum_i a_i\\,(x_i - \\hat{x}_i)
```

It is the weighted mean of the corrected clocks: the true phase ``x_i`` of
each clock less the filter's estimate ``\\hat{x}_i`` of it, `r.xhat`, with
the weights ``a_i`` of the scale
([Greenhall (2006), Eqs. (3.8) and (3.13)](@cite greenhall-2006-measurement-noise)).
The mean carries the reading noise of no clock, so its stability is that
of the scale itself. A measurement cannot give it: [`timescale_phase`](@ref)
reads the scale at one clock and adds the noise of that clock's reading.
With every `R = 0` the corrected clocks of the clocks read agree, and the
readout is the mean with its sign reversed.

# Arguments

- `r`: a scale from [`timescale`](@ref) with `algorithm = :greenhall`.
- `truth`: the true states from [`simulate_ensemble`](@ref);
  `truth[k, 1, i]` is the phase of clock ``i`` at sample ``k``.

# Examples

```julia
using SigmaTau, Random
ms = [ClockModel(; tau0 = 1.0, R = 1e-24, q1 = q, q2 = 1e-6q)
      for q in (1e-24, 4e-24, 1.6e-23)]
cc, truth = simulate_ensemble(ms, 20_000; reference = 0, rng = Xoshiro(1))
r = timescale(cc; algorithm = :greenhall, models = ms)
adev(implicit_mean(r, truth), [1])      # the scale itself
adev(timescale_phase(r), [1])           # read at clock 1, with its reading noise
```
"""
function implicit_mean(r::TimescaleResult, truth::AbstractArray{<:Real,3})
    r.algorithm === :greenhall ||
        throw(ArgumentError("implicit_mean: the implicit mean belongs to :greenhall, not :$(r.algorithm)"))
    size(truth, 1) == size(r.weights, 1) && size(truth, 3) == size(r.weights, 2) ||
        throw(ArgumentError("implicit_mean: truth is $(size(truth)), the scale has $(size(r.weights, 1)) samples and $(size(r.weights, 2)) clocks"))
    m = Float64[sum(r.weights[n, i] * (truth[n, 1, i] - r.xhat[n, i]) for i in axes(r.weights, 2)) for n in axes(r.weights, 1)]
    return PhaseData(m, r.comparisons.tau0; source="IEM")
end

# Median by sorting: Statistics.median is not trim-safe.
function _middle(v::Vector{Float64})
    s = sort(v)
    n = length(s)
    return isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1]) / 2
end

# Cap the weights at wmax and hand the excess to the clocks below the cap,
# in proportion to their own weight, until none is over.
function _cap_weights!(w::Vector{Float64}, wmax::Float64)
    capped = falses(length(w))
    while true
        over = w .> wmax
        any(over) || return w
        excess = sum(w[over]) - wmax * count(over)
        w[over] .= wmax
        capped .|= over
        free = .!capped .& (w .> 0)
        any(free) || return w
        w[free] .+= excess .* w[free] ./ sum(w[free])
    end
end

# Time variance ½E[(x₍ₙ₊₂ₘ₎ − 2x₍ₙ₊ₘ₎ + xₙ)²] = (mτ₀ σ_y(mτ₀))² of each clock (s²),
# from the second differences of the readings, by an N-cornered hat over the
# clocks and the reference. For two of them, half the mean square of the
# difference of their second differences is Sᵢⱼ = vᵢ + vⱼ, in which whatever
# they share cancels; the vᵢ are the least-squares solution weighted by
# n/Sᵢⱼ², the inverse variance of Sᵢⱼ. A clock the solution leaves nothing
# for takes half its smallest Sᵢⱼ. NaN where a clock has fewer than 8 triples
# in common with every other.
function _time_variances(r::AbstractMatrix{Float64}, m::Int)
    N, k = size(r)
    N > 2m || return fill(NaN, k)
    D = hcat(r[1+2m:N, :] .- 2 .* r[1+m:N-m, :] .+ r[1:N-2m, :], zeros(N - 2m))   # last: the reference
    S, A, b = fill(Inf, k + 1, k + 1), zeros(k + 1, k + 1), zeros(k + 1)
    for i in 1:k+1, j in i+1:k+1
        n, s = 0, 0.0
        for t in axes(D, 1)
            d = D[t, i] - D[t, j]
            isfinite(d) && (n += 1; s += d^2)
        end
        n >= 8 && s > 0 || continue         # s = 0: the reference is this clock
        S[i, j] = S[j, i] = s / 2n
        w = n / S[i, j]^2
        A[i, i] += w; A[j, j] += w; A[i, j] += w; A[j, i] += w
        b[i] += w * S[i, j]; b[j] += w * S[i, j]
    end
    use = [i for i in 1:k+1 if A[i, i] > 0]
    v = fill(NaN, k + 1)
    # The ridge keeps the solution where the pairs do not determine it.
    isempty(use) || (v[use] = (A[use, use] + 1e-12 * Diagonal(A[use, use])) \ b[use])
    for i in use
        v[i] > 0 || (v[i] = minimum(view(S, i, :)) / 2)
    end
    return v[1:k]
end

# Weiss & Weissert 1991 eq. 10: the exponential constant that best predicts a
# clock with white FM and random-walk FM whose ADEV has its minimum at τ_min.
_filter_constant(τmin::Real, τ₀::Real) = (-1 + sqrt(1 / 3 + 4 * (τmin / τ₀)^2 / 3)) / 2

# The stretch the starting values are taken from, and the usual spacing of
# the readings in it: the last of the first max(16, P/10) of the P samples
# with a reading, and the median number of samples from one to the next.
function _start(r::AbstractMatrix{Float64})
    rows = findall(n -> any(isfinite, view(r, n, :)), axes(r, 1))
    head = rows[1:min(length(rows), max(16, length(rows) ÷ 10))]
    return last(head), length(head) < 2 ? 1 : round(Int, _middle(Float64.(diff(head))))
end

# Each clock's own τ_min: the lag where its individual ADEV is smallest among
# the octaves of the spacing `c` up to a tenth of the record. Beyond them few
# differences decide, and an octave counts only while the ADEV still falls
# from one octave to the next: a clock whose minimum the record does not
# reach keeps the longest lag, and no other clock takes one by chance.
function _estimate_tau_min(cc::ClockComparisons, c::Int)
    N, k = size(cc.x)
    best, tmin, falling = fill(Inf, k), fill(c * cc.tau0, k), falses(k)
    m = c
    while N - 2m >= 16c
        v = _time_variances(cc.x, m)
        for i in 1:k
            s = sqrt(v[i]) / m
            falling[i] = isfinite(s) && s < best[i] && (10m <= N || falling[i])
            falling[i] && (best[i] = s; tmin[i] = m * cc.tau0)
        end
        m *= 2
    end
    return tmin
end

# Starting prediction-error variances ε_i²(τ) = τ² σ_y,i²(τ) (eq. 7) at the
# spacing τ = cτ₀, each clock's from the `stop` samples that begin at its
# first reading; a clock with no usable estimate takes the largest of the others.
function _initial_variance(cc::ClockComparisons, stop::Int, c::Int)
    N, k = size(cc.x)
    v = fill(NaN, k)
    for i in 1:k
        a = findfirst(isfinite, view(cc.x, :, i))
        a === nothing || (v[i] = _time_variances(view(cc.x, a:min(N, a + stop - 1), :), c)[i])
    end
    ok = filter(s -> isfinite(s) && s > 0, v)
    isempty(ok) && return ones(k)
    return Float64[isfinite(s) && s > 0 ? s : maximum(ok) for s in v]
end

# Least-squares slope of clock i's readings over the `stop` samples from
# sample `a`; NaN with fewer than two readings there.
function _slope(cc::ClockComparisons, i::Int, a::Int, stop::Int)
    rows = a:min(size(cc.x, 1), a + stop - 1)
    ok = isfinite.(cc.x[rows, i])
    count(ok) >= 2 || return NaN
    t = (rows[ok] .- 1) .* cc.tau0
    return sum((t .- sum(t) / length(t)) .* cc.x[rows[ok], i]) / sum(abs2, t .- sum(t) / length(t))
end

# Starting frequencies against the scale of the clocks read at the first
# sample: each one's slope over the first `stop` samples less the mean of the
# slopes with the weights `w`, so clocks with frequency offsets start tracking
# them (the weighted mean of the offsets is the scale's own). Zero for a clock
# without a slope.
function _initial_frequency(cc::ClockComparisons, stop::Int, w::Vector{Float64})
    s = Float64[isfinite(cc.x[1, i]) ? _slope(cc, i, 1, stop) : NaN for i in axes(cc.x, 2)]
    fit = isfinite.(s)
    return Float64[fit[i] ? s[i] - sum(w[fit] .* s[fit]) / sum(w[fit]) : 0.0 for i in eachindex(s)]
end

# Starting frequency against the scale of clock i, which joins at sample n
# for the first time: its slope over the `stop` samples from n less the
# scale's against the reference
# there, which is the weighted mean, over the clocks in the average, of each
# one's slope less its frequency against the scale. The reference cancels.
function _joining_frequency(cc::ClockComparisons, i::Int, n::Int, stop::Int, w::Vector{Float64}, y::Vector{Float64})
    s, sw = 0.0, 0.0
    for j in eachindex(w)
        sⱼ = w[j] > 0 && j != i ? _slope(cc, j, n, stop) : NaN
        isfinite(sⱼ) && ((s, sw) = (s + w[j] * (sⱼ - y[j]), sw + w[j]))
    end
    sᵢ = _slope(cc, i, n, stop)
    return isfinite(sᵢ) && sw > 0 ? sᵢ - s / sw : 0.0
end

# Steady-state frequency variance P of the AT2 filter (eq 16 for any τ): with
# s = σ_β²τ, the fixed point of P̂ = P + s, P = σα² P̂/(σα² + P̂), i.e.
# P² + sP − σα²s = 0. Units: σα² and s are fractional frequency squared.
_at2_steady(σα²::Float64, s::Float64) = (-s + sqrt(s^2 + 4 * σα² * s)) / 2

"""
    timescale(cc::ClockComparisons; algorithm=:at1, models=nothing, pivot=1,
              tau_min=nothing, memory=1.728e6, bias=:paper, wmax=1.0,
              outlier=Inf, steps=false, step_threshold=5.0, lmax=nothing,
              step_inflate=4.0) → TimescaleResult

A time scale from the clocks of `cc`: the phase and frequency of each
clock against the scale, and its weight, at every sample.

| `algorithm`    | Scale                                     | Weights                       |
|:---------------|:------------------------------------------|:------------------------------|
| `:ewfa`        | equal-weights frequency averaging         | equal                         |
| `:at1`         | NIST AT1                                  | from the prediction errors    |
| `:at2`         | AT1 with a frequency variance             | as `:at1`                     |
| `:greenhall`   | Kalman filter with covariance reduction   | of the implicit ensemble mean |
| `:jones_tryon` | the same filter without the reduction     | not reported                  |
| `:stein`       | a Kalman filter for each clock difference | inverse noise variances       |

The first three learn what they need from the comparisons. The last three
take one [`ClockModel`](@ref) per clock in `models`, which
[`clock_models`](@ref) fits from `cc`. [Time scales](@ref) compares the six
in simulation.

# Averaging scales

With ``x_i`` the phase of clock ``i`` against the scale, ``y_i`` its
frequency and ``z_i`` its reading against the reference, each clock
predicts its next phase, and the scale is placed where the weighted
prediction errors sum to zero
([Weiss and Weissert (1991), Eqs. (1)–(2)](@cite weiss-weissert-1991-at2),
without the drift term of Eq. (1)):

```math
\\hat{x}_i = x_i(t) + y_i(t)\\,τ, \\qquad
x_i(t + τ) = z_i(t + τ) + \\sum_j w_j \\bigl[\\hat{x}_j - z_j(t + τ)\\bigr] .
```

A clock without a reading keeps its prediction and takes no part in the
sum. A clock joins at its first reading, and again at the first reading
after a gap: its phase is then taken from the scale, and it takes part in
the sum from the next sample, so the scale does not step when a clock
returns. Where every clock read at a sample returns from a gap, they all
take part. The sources give no procedure for a clock that returns, and
this one is the package's own.
[Weiss and Weissert (1991), §1](@cite weiss-weissert-1991-at2) name the
aim: clocks are added to AT1 and removed from it with minimal impact on
the scale.

In these three scales the interval ``τ`` of the equations is the time since
the last sample with a reading, as a clock's frequency comes from its
latest two readings
([Weiss and Weissert (1991), §2.8](@cite weiss-weissert-1991-at2)): clocks
read together at every tenth sample give the scale of those readings alone,
and the constants below are taken at that spacing. The clocks have to be
read at the same samples; one read less often than the others joins anew at
each of its readings and never has a weight.

The scale starts at the mean of the first readings. Each clock starts from
its own first readings, as many samples as a tenth of the record's
readings. The frequency of a clock read at the first sample starts at the
slope of a line through them, less the mean slope of those clocks under
their starting weights. A clock that joins later starts at its slope less that of the scale over the same
samples, which is the weighted mean of the other clocks' slopes less their
frequencies against the scale, so the reference of the readings cancels.
The rate of the scale then changes by the error of that frequency times
the weight of the clock.

`:ewfa` gives every clock the same weight and takes ``y_i`` as the clock's
last frequency step, ``[x_i(t) - x_i(t - τ)]/τ``
([Kettering (2023), §3.1](@cite kettering-2023-clock-ensemble)). With every
clock read at every sample the scale is the mean of the clocks' phases. The
thesis does not treat missing readings. Under the rule above the scale
keeps its phase when a clock drops out, and its frequency shifts by the
noise of one frequency step until the clock returns.

`:at1` is the NIST scale as
[Weiss and Weissert (1991), Eqs. (3)–(10)](@cite weiss-weissert-1991-at2)
give it. The weights are the inverse variances of the prediction errors,

```math
w_i = \\frac{ε_x^2}{ε_i^2}, \\qquad ε_x^2 = \\Bigl(\\sum_i ε_i^{-2}\\Bigr)^{-1} ,
```

where ``ε_i^2`` is an exponential average, of time constant `memory`, of
the square of

```math
\\hat{ε}_i = |\\hat{x}_i - x_i| + K_i, \\qquad K_i = \\frac{2 ε_x^2}{\\sqrt{2π ε_i^2}} .
```

``K_i`` makes up for the clock being part of the scale it is compared
with. ``ε_i^2`` starts at ``τ^2 σ_{y,i}^2(τ)`` (§2.7 of the paper),
estimated from those first readings of the clock, so a clock that joins
late has its weight at once. The estimate is an N-cornered hat over the
clocks and the reference
([Gray and Allan (1974)](@cite gray-1974-three-cornered-hat);
[SP1065 §10.14](@cite riley-2008-sp1065)), so that the noise of the
reference and of the other clocks stays out of it; each pair has the weight
of the inverse variance of its sum, which is the package's own. The
frequency is an exponential average of the frequency steps,

```math
y_i(t + τ) = \\frac{1}{m_i + 1} \\Bigl[\\frac{x_i(t + τ) - x_i(t)}{τ} + m_i\\,y_i(t)\\Bigr],
\\qquad m_i = \\frac{1}{2} \\Bigl[\\Bigl(\\frac{1}{3} + \\frac{4 τ_{\\min,i}^2}{3 τ^2}\\Bigr)^{1/2} - 1\\Bigr] ,
```

with ``τ_{\\min,i}`` the averaging time at which the Allan deviation of the
clock is smallest.

`:at2` keeps the weights of AT1 and estimates each frequency with a Kalman
filter on the frequency steps ``\\bar{y}_i``
([Weiss and Weissert (1991), Eqs. (11)–(18)](@cite weiss-weissert-1991-at2)).
The white FM of a step has the variance ``σ_α^2 = ε_i^2/τ^2``, and the
random walk of the frequency gains ``σ_β^2 τ`` per interval:

```math
\\hat{P} = P + σ_β^2 τ, \\qquad
y_i \\leftarrow \\frac{σ_α^2\\,y_i + \\hat{P}\\,\\bar{y}_i}{σ_α^2 + \\hat{P}}, \\qquad
P \\leftarrow \\frac{σ_α^2\\,\\hat{P}}{σ_α^2 + \\hat{P}} .
```

``σ_β^2 τ`` is fixed for each clock by ``m_i(m_i + 1) = σ_α^2/(σ_β^2 τ)``
at the start, where ``P`` takes its steady state. As ``σ_α^2`` follows
``ε_i^2``, the memory of the filter adapts to the clock; in steady state it
is the filter of AT1.

# Frequency steps

With `steps = :best` or `:real`, `:at2` looks for a step in the frequency
of each clock after every sample, after
[Weiss and Weissert (1991), §4](@cite weiss-weissert-1991-at2). For windows
of ``s = 2`` to `lmax` − 1 intervals that end at the current sample, the
mean frequency of the clock over the window is compared with its estimate
``y_i`` at the start of the window. The variance of the difference is
taken as

```math
P_{-s} + \\frac{σ_{α,i}^2 + σ_{αx}^2}{s}
+ \\bigl(σ_{β,i}^2 τ + σ_{βx}^2 τ\\bigr) \\frac{(s + 1)(2s + 1)}{6s} ,
```

the variance of the estimate at the start of the window, the white FM of a
mean over ``s`` intervals and the random walk of the frequency across
them, the last two for the clock and for the scale, whose levels are
``σ_{αx}^2 = (\\sum_i σ_{α,i}^{-2})^{-1}`` and likewise ``σ_{βx}^2``
(Eqs. (21)–(22)). The largest difference beyond `step_threshold` standard
deviations marks a step at the start of its window. No window starts
before the sample at which the clock last joined the scale, or returned to
it after a step.

The scale is then run again from the step with the weight of that clock
at zero and its ``P`` raised to ``σ_α^2``, and no further step is looked
for until the rerun has reached the current sample. The clock stays out for
``τ_{\\min,i}``. Meanwhile it is still read against the scale and its
frequency and ``ε_i^2`` go on being estimated, so that it returns with its
new frequency, and ``ε_i^2`` is multiplied by `step_inflate` at its return.
`:best` replaces the scale from the step on, which suits a recorded
history. `:real` keeps every value already issued and corrects only what
follows, as a scale that runs in real time must.

The test differs from the paper's in four ways. The variance above stands
in place of that of Eq. (20). The windows end at the current sample and are
at least two intervals long, where the average frequency of Eq. (19) ends
one measurement earlier. The default threshold is 5 standard deviations,
where the paper takes 4. And the paper gives neither the amount by which
``P`` is raised nor that by which ``ε_i^2`` is. [Frequency-step detection](@ref)
gives the steps caught and the false alarms in simulation.

# Kalman scales

`:greenhall` runs one Kalman filter over the states of all the clocks, with
the transition and process noise of each model and the differences of the
readings as measurements. After each measurement update the covariance
``P`` is reduced in its phase part ``P_{xx}``:

```math
P \\leftarrow S P S^\\mathsf{T}, \\qquad
S = \\begin{pmatrix} I - \\mathbf{1} a^\\mathsf{T} & 0 \\\\ 0 & I \\end{pmatrix}, \\qquad
a = \\frac{P_{xx}^{-1} \\mathbf{1}}{\\mathbf{1}^\\mathsf{T} P_{xx}^{-1} \\mathbf{1}}
```

([Greenhall (2006), Eqs. (3.11)–(3.12) and (5.1)](@cite greenhall-2006-measurement-noise)).
The reduction leaves the frequency estimates as they were, and the weights
``a`` are those that minimize the variance of each step of the scale
([Greenhall (2003), Theorems 1 and 2](@cite greenhall-2003-jones-tryon);
Theorem B.2 of the 2006 paper with noise on the readings). With every
`R = 0` this is the reduced Kalman scale of the 2003 paper, and with
`R > 0` the IEMx algorithm of the 2006 paper. The weights returned are
``a``, those of the implicit ensemble mean of [`implicit_mean`](@ref).

A clock without a reading takes no part in the measurement or in ``a``,
which is then taken over the clocks read, the procedure of
[Greenhall (2006), §10](@cite greenhall-2006-measurement-noise). The clock
runs on its prediction, and the scale does not step when its readings
return.

A clock first read after the start begins as the clock of most weight read
with it, plus a line through its first 16 readings against that clock. The
covariance of the line is widened, so that the readings to come set the
phase and frequency of the clock. The scale does not step there, and the
reference of the readings takes no part in the start. The sources give no
procedure, and this one is the package's own:
[Greenhall (2003), §5](@cite greenhall-2003-jones-tryon) names the
insertion of a new clock without a phase or frequency step as a practical
difficulty of a Kalman scale.

Noise on the readings leaves the corrected clocks ``x_i - \\hat{x}_i``
unequal and their mean out of reach of a measurement (§4 of the 2006
paper). The scale is therefore read at clock `pivot`, as the reading of
that clock less its estimated phase, and it carries the white PM of that
reading. With every `R = 0` all pivots give the same scale.

The filter starts each clock read at the first sample from a line through
its first 16 readings. An error in those starting frequencies stays in the
scale as a constant frequency offset
([Greenhall (2007), §2](@cite greenhall-2007-measurement-noise-corrections)).
The drift of a model with drift starts at zero, as in the initialization of
[Jones and Tryon (1983), Fig. 1](@cite jones-tryon-1983-estimating-time),
with the variance that a parabola through those readings gives it. The
filter observes the differences between the drifts of the clocks and not
their common part, so a drift fitted to the first readings would stay in
the scale as a constant drift; started at zero, the drift of the scale is a
weighted mean of the drifts of the clocks.
Models with flicker FM bring their Markov states into the filter, which at
`R = 0` is the algorithm of
[Davis et al. (2005), Eqs. (3.1)–(3.6) and (4.1)](@cite davis-greenhall-stacey-2005-flicker);
the frequency returned for a clock includes them. Markov states for
flicker PM belong to the readings and stay out of the phases that the
scale is formed from. The papers take the noise of a reading as white, so
that case is this package's own; see [Flicker states](@ref).

`:jones_tryon` is the same filter without the reduction. Its scale is the
one that [Jones and Tryon (1983), §4](@cite jones-tryon-1983-estimating-time)
read from their filter, the reading of a clock less its estimated phase,
which [Greenhall (2003), Eq. (10)](@cite greenhall-2003-jones-tryon) calls
the raw Kalman scale and finds unstable at short averaging times. The
frequency estimates are those of `:greenhall`. The weights are returned as
`NaN`: those implied by the Kalman gain (Eq. (19) of the 2003 paper) are
not worked out. Against the 1983 paper, the process noise is the
integrated one of [`process_noise`](@ref) where the paper's is diagonal
(its Eq. (2.2)), the differences share the reading noise of the clock they
are taken against where the paper's errors are independent (Eq. (2.6)),
and the paper's detection and correction of clock errors
(Eqs. (4.6)–(4.13)) is left out. A missing reading leaves its difference
out of the measurement, as in §2 of the paper, which gives no rule for the
scale. When the readings of a clock return after a gap, the scale takes a
share of the clock's prediction error and steps. A clock first read after
the start begins as in `:greenhall`.

`:stein` is the scale of [Stein (2003)](@cite stein-2003-demystified), the
KAS-2 algorithm that
[Levine (2012), §XVI](@cite levine-2012-statistical-modeling) describes.
The phase, the frequency and, where a model has `q3 > 0`, the aging of
each clock take independent shocks (Eq. (3.1)), here of variance ``q_1 τ``,
``q_2 τ`` and ``q_3 τ``. Each clock is differenced against clock `pivot`,
and each difference has its own Kalman filter (Eq. (5.8)). The shocks of
the pivot follow from the constraint that the weighted shocks of all the
clocks sum to zero, with one set of weights for each state
(Eqs. (4.1), (6.2) and (7.1)–(7.2)), and the other clocks are the pivot
plus their differences. The weights are inversely proportional to the
shock variances
([Erickson et al. (2011), §3.1](@cite erickson-2011-symmetricom-time-scale))
and capped at `wmax`; those returned are the phase weights. The scale is
read at the pivot as for `:greenhall`. The paper runs a filter for every
pair of clocks; with the filters of one pivot only, the scale depends on
the pivot. Models with flicker states are refused; `clock_model` makes
one without them with `flicker = false`. A flicker level too low for
states counts as white noise, as in the other filters.

The paper starts from given values of the states (after Eq. (7.3)). Here
the difference of each clock starts from the estimate of its filter at the
first reading of the clock, smoothed back over the clock's first tenth of
the record, 16 readings at least. The scale starts as the weighted mean of
the clocks read at the first sample, in phase, frequency and aging, each
with the weights of that state. The procedure is the package's own.

The paper does not treat missing readings. Here a clock without a reading
runs on its prediction, and at its first reading after a gap it takes no
part in the shocks of the pivot, because its innovation spans the gap; the
scale does not step. While the pivot itself is unread no difference is
measured, every clock runs on its prediction, and the scale steps when the
pivot returns.

!!! note "What the learned weights are optimal for"
    `:at1` and `:at2` weigh a clock by its prediction error over one
    sample interval. At averaging times far beyond it, among clocks whose
    long-term noise differs in kind, the scale can be less stable than its
    best clock. With two clocks the comparisons cannot tell which is the
    better, and the weights stay equal when one is read against the other.

!!! note "Missing readings"
    No scale steps when a clock drops out. When the readings of a clock
    return, `:jones_tryon` steps, and so does `:stein` when the clock is
    its pivot: give `:stein` a pivot that is read throughout. The other
    scales do not step. [Clocks that drop out and return](@ref) gives the
    sizes in simulation.

# Arguments

- `cc`: the clocks, read against one reference. The scale starts at the
  first sample at which a clock is read; the samples before it are `NaN`
  in the result.

# Keywords

- `algorithm`: one of the six scales of the table.
- `models`: one [`ClockModel`](@ref) per clock, with the `tau0` of `cc`,
  for the Kalman scales. At most one may be free of noise.
- `pivot`: the clock at which a Kalman scale is read. Where it has no
  reading, the first clock read stands in. `:stein` needs its pivot read at
  the scale's first sample, and steps when the pivot returns from a gap.
- `tau_min`: ``τ_{\\min}`` in seconds for `:at1` and `:at2`, one number or
  one per clock. By default, the τ of the smallest Allan deviation of each
  clock among the octaves of the interval between readings, estimated from
  the record as the starting ``ε_i^2`` is. Beyond a tenth of the record an
  octave counts only while the deviation still falls.
- `memory`: time constant in seconds of the average behind ``ε_i^2``. The
  default is 20 days, the paper's value for cesium clocks; a quarter of the
  record is used where that is shorter.
- `bias`: `:paper` for ``K_i`` above, or `:exact` to divide the squared
  prediction error by ``1 - w_i`` instead, the correction of
  [Levine (2012), Eq. (27)](@cite levine-2012-statistical-modeling), which
  is exact for independent clocks at the weights above.
- `wmax`: the largest weight a clock may take in `:at1`, `:at2` and
  `:stein`; the excess goes to the other clocks in proportion to their
  weights. The NIST scale limits a weight to 0.3 (§IX and Eq. (27) of
  that review).
- `outlier`: in `:at1` and `:at2`, leave a clock out of a sample when its
  prediction error is more than `outlier` times ``ε_i`` from the median
  error of the clocks read, unless that leaves none. Off by default.
- `steps`: `false`, `:best` (also `true`) or `:real`, for `:at2`.
- `step_threshold`: standard deviations that a frequency difference must
  pass to count as a step.
- `lmax`: one more than the longest window of the step test, in intervals;
  by default ``τ_{\\min,i}/τ`` for each clock, and at least 4.
- `step_inflate`: factor on ``ε_i^2`` when a clock returns after a step;
  4 doubles ``σ_α``.

# Examples

```jldoctest
julia> using SigmaTau

julia> cc = ClockComparisons([0.0 2.0; 1.0 5.0; 2.0 8.0; 3.0 11.0] .* 1e-9, 1.0);

julia> r = timescale(cc; algorithm = :ewfa);

julia> r.x[:, 1] ≈ -r.x[:, 2]      # the scale is midway between two clocks
true
```

Three simulated clocks under the reduced Kalman scale, which weighs them
by the inverse of their white FM:

```julia
using SigmaTau, Random
ms = [ClockModel(; tau0 = 1.0, q1 = q, q2 = 1e-6q) for q in (1e-24, 4e-24, 1.6e-23)]
cc, _ = simulate_ensemble(ms, 20_000; reference = 0, rng = Xoshiro(1))
r = timescale(cc; algorithm = :greenhall, models = ms)
r.weights[end, :]              # 0.762, 0.190, 0.048
adev(timescale_phase(r))       # the scale against the ideal reference
```
"""
function timescale(cc::ClockComparisons; algorithm::Symbol=:at1,
                   tau_min::Union{Nothing,Real,AbstractVector{<:Real}}=nothing,
                   memory::Real=20 * 86400.0, bias::Symbol=:paper, wmax::Real=1.0,
                   outlier::Real=Inf, models::Union{Nothing,AbstractVector{ClockModel}}=nothing,
                   pivot::Integer=1, steps::Union{Bool,Symbol}=false, step_threshold::Real=5.0,
                   lmax::Union{Nothing,Integer}=nothing, step_inflate::Real=4.0)
    algorithm in (:ewfa, :at1, :at2, :greenhall, :jones_tryon, :stein) ||
        throw(ArgumentError("timescale: algorithm must be :ewfa, :at1, :at2, :greenhall, :jones_tryon or :stein, got :$algorithm"))
    bias in (:paper, :exact) || throw(ArgumentError("timescale: bias must be :paper or :exact, got :$bias"))
    mode = steps === true ? :best : steps === false ? :none : steps
    mode in (:none, :best, :real) ||
        throw(ArgumentError("timescale: steps must be false, :best or :real, got :$steps"))
    mode !== :none && algorithm !== :at2 && throw(ArgumentError("timescale: steps needs algorithm :at2"))
    kalman = algorithm in (:stein, :greenhall, :jones_tryon)
    kalman || models === nothing ||
        throw(ArgumentError("timescale: :$algorithm takes no models; they are for :greenhall, :jones_tryon and :stein"))
    tau_min === nothing || all(>(0), tau_min) || throw(ArgumentError("timescale: tau_min must be positive, got $tau_min"))
    memory >= 0 || throw(ArgumentError("timescale: memory must be ≥ 0, got $memory"))
    step_threshold > 0 || throw(ArgumentError("timescale: step_threshold must be positive, got $step_threshold"))
    N, k = size(cc.x)
    N >= 3 || throw(ArgumentError("timescale: needs three or more samples"))
    n₁ = findfirst(n -> any(isfinite, view(cc.x, n, :)), 1:N)
    n₁ === nothing && throw(ArgumentError("timescale: no clock has a reading"))
    if n₁ > 1       # the scale starts at the first sample with a reading
        r = timescale(ClockComparisons(cc.x[n₁:N, :], cc.tau0; names=cc.names, reference=cc.reference);
                      algorithm, tau_min, memory, bias, wmax, outlier, models, pivot, steps, step_threshold, lmax, step_inflate)
        lead(A, v) = vcat(fill(v, n₁ - 1, k), A)
        return TimescaleResult(algorithm, cc, collect((0:N-1) .* cc.tau0), lead(r.x, NaN), lead(r.y, NaN),
                               lead(r.weights, algorithm === :jones_tryon ? NaN : 0.0), lead(r.error, NaN), lead(r.xhat, NaN))
    end
    if kalman
        models === nothing && throw(ArgumentError("timescale: :$algorithm needs models, one ClockModel per clock"))
        1 <= pivot <= k || throw(ArgumentError("timescale: pivot $pivot is not one of the $k clocks"))
        count(_noiseless, models) <= 1 || throw(ArgumentError(
            "timescale: two or more clock models have no noise, so the filter cannot weigh their readings"))
        algorithm === :stein && return _stein(cc, models, Int(pivot), Float64(wmax))
        return _reduced_kalman(cc, models, Int(pivot), algorithm)
    end
    k * wmax >= 1 || throw(ArgumentError("timescale: wmax below 1/$k leaves weights that cannot sum to 1"))
    at2 = algorithm === :at2
    at1 = algorithm === :at1 || at2
    τ = cc.tau0
    r = cc.x
    stop, c = _start(r)         # c: the usual spacing of the readings, in samples
    τmin = tau_min === nothing ? (at1 ? _estimate_tau_min(cc, c) : fill(τ, k)) :
           tau_min isa Real ? fill(Float64(tau_min), k) : Float64.(tau_min)
    length(τmin) == k || throw(ArgumentError("timescale: $(length(τmin)) tau_min values for $k clocks"))
    mf = at1 ? Float64[_filter_constant(t, c * τ) for t in τmin] : zeros(k)
    Nτ = max(1.0, min(memory, N * τ / 4) / (c * τ))
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    ε² = at1 ? _initial_variance(cc, stop, c) : ones(k)
    # Start: the scale is the plain mean of the first readings.
    ok = isfinite.(r[1, :])
    u = -sum(r[1, ok]) / count(ok)
    X[1, :] = r[1, :] .+ u
    # The interval of eqs 7-15 is the time since the last sample with a reading,
    # `prev`; a clock joins at a reading that follows none of its own there.
    seen, prev = copy(ok), 1
    w = ok ./ ε²
    y = _initial_frequency(cc, stop, at1 ? _cap_weights!(w ./ sum(w), max(Float64(wmax), 1 / count(ok))) : w)
    # AT2: the random-walk FM impulse σ_β²τ of each clock from its τ_min through
    # the steady state m(m+1) = σα²/(σ_β²τ) (eq 18), and the frequency variance P
    # at its steady state (eq 16), both from the starting σα² = ε²/τ² at τ = cτ₀.
    σβτ = ε² ./ (c * τ)^2 ./ (mf .* (mf .+ 1))
    Pf = Float64[_at2_steady(ε²[i] / (c * τ)^2, σβτ[i]) for i in 1:k]
    W[1, :] = ok ./ count(ok)
    # Step detector state (AT2, eqs 19-22): the state after each sample, kept
    # to rerun the scale from a detected step; `blocked[i]` is the last sample
    # clock i is kept out of the scale; no detection until `replay`.
    Lmax = [lmax === nothing ? max(4, round(Int, t / (c * τ))) : Int(lmax) for t in τmin]
    ε²h, yh, Pfh, prevh = zeros(N, k), zeros(N, k), zeros(N, k), ones(Int, N)
    blocked, replay = zeros(Int, k), 0
    joined = ones(Int, k)       # the sample each clock last joined at: no step window reaches behind it
    # :real keeps what was issued: rows are copied out once, after a rerun
    # continued from the last issued phase (the correction steers forward).
    real = mode === :real
    Xo, Yo, Wo, Eo = real ? (zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)) : (X, Y, W, E)
    issued, pending, shift = 0, false, 0.0
    n = 2
    while n <= N
        if real && n - 1 > issued
            if pending     # phase continuity with the last readings issued
                p = prevh[issued]
                j = findall(isfinite, view(r, p, :))
                shift = isempty(j) ? 0.0 : sum(Xo[p, j] .- X[p, j]) / length(j)
                pending = false
            end
            Xo[n-1, :], Yo[n-1, :], Wo[n-1, :], Eo[n-1, :] = X[n-1, :] .+ shift, Y[n-1, :], W[n-1, :], E[n-1, :]
            issued = n - 1
        end
        ε²h[n-1, :], yh[n-1, :], Pfh[n-1, :], prevh[n-1] = ε², y, Pf, prev
        for i in 1:k
            n == blocked[i] + 1 && (ε²[i] *= step_inflate)   # back in, with σα raised
        end
        x̂ = X[n-1, :] .+ y .* τ
        read = isfinite.(r[n, :])
        use = read .& seen .& (n .> blocked)
        # Every clock read is back from a gap. One that has never been read
        # with the scale has no place in it yet (its x̂ is NaN).
        had = read .& isfinite.(x̂)
        any(use) || (use = had .& (n .> blocked))
        if at1 && isfinite(outlier) && n > 2 && any(use)
            # The median prediction error places the scale for the test, so
            # that one clock's jump cannot carry it.
            d = x̂ .- r[n, :]
            keep = use .& (abs.(d .- _middle(d[use])) .<= outlier .* sqrt.(ε²))
            any(keep) && (use = keep)
        end
        if !any(use)
            X[n, :], Y[n, :] = x̂, y
            any(had) && ((seen, prev) = (had, n))
            n += 1
            continue
        end
        wu = at1 ? use ./ ε² : Float64.(use)
        wu ./= sum(wu)
        εx² = at1 ? 1 / sum(1 ./ ε²[use]) : 0.0
        _cap_weights!(wu, max(Float64(wmax), 1 / count(use)))
        u = sum(wu[use] .* (x̂[use] .- r[n, use]))
        for i in 1:k
            if read[i] && !seen[i] && !use[i]
                X[n, i] = r[n, i] + u
                joined[i] = n
                isfinite(x̂[i]) || (y[i] = _joining_frequency(cc, i, n, stop, wu, Y[n-1, :]))
            elseif read[i]
                X[n, i] = r[n, i] + u
                E[n, i] = X[n, i] - x̂[i]
                # A clock kept out after a step is still estimated, so that it
                # comes back with its new frequency; its bias K = wδ is zero.
                # One the outlier test left out keeps its variance and frequency.
                est = use[i] || n <= blocked[i]
                if at1 && est
                    ε̂² = bias === :paper ? (abs(E[n, i]) + (use[i] ? 2εx² / sqrt(2π * ε²[i]) : 0.0))^2 :
                         E[n, i]^2 / (1 - min(wu[i], 0.9))
                    ε²[i] = (ε̂² + Nτ * ε²[i]) / (1 + Nτ)
                end
                if est
                    m = n - prev
                    ŷ = (X[n, i] - X[prev, i]) / (m * τ)
                    if at2          # Kalman update of y (eqs 13-15), σα² = ε²/τ²
                        σα² = ε²[i] / (c * m * τ^2)
                        P̂ = Pf[i] + σβτ[i] * m / c
                        y[i] = (σα² * y[i] + P̂ * ŷ) / (σα² + P̂)
                        Pf[i] = σα² * P̂ / (σα² + P̂)
                    else
                        mᵢ = at1 ? _filter_constant(τmin[i], m * τ) : 0.0
                        y[i] = (mᵢ * y[i] + ŷ) / (mᵢ + 1)
                    end
                end
            else
                X[n, i] = x̂[i]
            end
        end
        seen, prev = read, n
        Y[n, :] = y
        W[n, :] = wu
        if mode !== :none && n > replay
            n0, i0 = _at2_step(X, yh, Pfh, ε², σβτ, use, max.(blocked .+ 1, joined), Lmax, n, c, c * τ, step_threshold)
            if i0 > 0       # a frequency step: rerun from its start without that clock
                blocked[i0] = n0 + max(1, round(Int, τmin[i0] / τ))
                replay = n
                pending = real
                ε², y, Pf = ε²h[n0, :], yh[n0, :], Pfh[n0, :]
                Pf[i0] = ε²[i0] / (c * τ)^2    # P raised at the step
                prev = prevh[n0]
                seen = isfinite.(r[prev, :])
                joined .= min.(joined, n0)
                n = n0 + 1
                continue
            end
        end
        n += 1
    end
    if real && N > issued       # the last row, when no rerun is pending for it
        if pending
            p = prevh[issued]
            j = findall(isfinite, view(r, p, :))
            shift = isempty(j) ? 0.0 : sum(Xo[p, j] .- X[p, j]) / length(j)
        end
        for m in issued+1:N
            Xo[m, :], Yo[m, :], Wo[m, :], Eo[m, :] = X[m, :] .+ shift, Y[m, :], W[m, :], E[m, :]
        end
    end
    Yo[isnan.(Xo)] .= NaN       # no frequency before a clock's first reading
    return TimescaleResult(algorithm, cc, collect((0:N-1) .* τ), Xo, Yo, Wo, Eo, Xo)
end

# AT2 frequency-step test (Weiss & Weissert 1991, eqs 19-22), run on the clocks
# kept at sample n: for each window of s = 2 … Lmax−1 intervals (τ, `c` samples
# each), the clock's mean frequency against the scale over the window, y_avg = (x[n] − x[n−sc])/(sτ),
# is compared with the frequency estimate y₋ₗ before it. In place of the
# variance of the paper's eq 20, that of the difference, in fractional frequency²:
#   P₋ₗ + (σα,i² + σαx²)/s + (σβτ,i + σβτx)(s+1)(2s+1)/(6s)
# the estimate's variance, the white-FM noise of the mean over s intervals (the
# clock's and the scale's), and the random walk of the frequency over the window
# (the mean of a random walk of impulse σβτ per interval over s intervals has
# variance σβτ (s+1)(2s+1)/(6s)). σαx² and σβτx combine the clocks' levels as in eqs 21-22,
# with σβτ the paper's σβ²τ. Returns (start sample, clock) of the window with
# the largest ratio above `thr`, or (0, 0).
function _at2_step(X, yh, Pfh, ε², σβτ, use, from, Lmax, n, c, τ, thr)
    σαx⁻, σβx⁻ = 0.0, 0.0
    for j in eachindex(use)
        use[j] || continue
        σαx⁻ += τ^2 / ε²[j]
        σβx⁻ += 1 / σβτ[j]
    end
    best, n0b, ib = thr, 0, 0
    for i in eachindex(use)
        use[i] || continue
        σα² = ε²[i] / τ^2
        for s in 2:min(Lmax[i] - 1, (n - from[i]) ÷ c)
            n0 = n - s * c
            var = Pfh[n0, i] + (σα² + 1 / σαx⁻) / s + (σβτ[i] + 1 / σβx⁻) * (s + 1) * (2s + 1) / (6s)
            z = abs((X[n, i] - X[n0, i]) / (s * τ) - yh[n0, i]) / sqrt(var)
            z > best && ((best, n0b, ib) = (z, n0, i))
        end
    end
    return n0b, ib
end

# IEM weights a (Σa = 1) minimizing aᵀ Pxx a: Pxx⁻¹1 / 1ᵀPxx⁻¹1. Clocks with
# noiseless readings (`quiet`) differ by known offsets, so only their total
# weight is determined: it is found with one of them standing for the group,
# and split among them as `split` (the gain's). `nothing` when Pxx is singular
# for another reason.
function _iem_weights(Pxx::Matrix{Float64}, quiet::AbstractVector{Bool}, split::Vector{Float64})
    g = findall(quiet)
    keep = [i for i in axes(Pxx, 1) if !(i in g) || i == first(g)]
    length(g) < 2 && (keep = collect(axes(Pxx, 1)))
    P = Pxx[keep, keep]
    s = maximum(diag(P))
    s > 0 || return nothing
    C = cholesky(Symmetric(P ./ s); check=false)
    (issuccess(C) && minimum(diag(C.U)) > 1e-6) || return nothing
    w = C \ ones(length(keep))
    a = zeros(size(Pxx, 1))
    a[keep] = w ./ sum(w)
    if length(g) >= 2
        t = sum(split[g])
        a[g] = a[first(g)] .* (t > 0 ? split[g] ./ t : fill(1 / length(g), length(g)))
    end
    return a
end

function _reduced_kalman(cc::ClockComparisons, models::AbstractVector{ClockModel}, pivot::Int, algorithm::Symbol)
    N, k = size(cc.x)
    length(models) == k || throw(ArgumentError("timescale: $(length(models)) models for $k clocks"))
    τ = cc.tau0
    all(m -> m.tau0 ≈ τ, models) || throw(ArgumentError("timescale: the models' tau0 must be the comparisons' τ₀"))
    ns = maximum(_nstates, models)
    # State: the k phases, then the k frequencies, then each clock's further
    # states (drift, Markov states of the flicker) in its own order.
    idx(i, s) = (s - 1) * k + i
    Φ = Matrix{Float64}(I, ns * k, ns * k)
    Q = zeros(ns * k, ns * k)
    X̂ = zeros(ns * k)
    P = zeros(ns * k, ns * k)
    G, F = zeros(k, ns * k), zeros(k, ns * k)     # each clock's reading and frequency in state terms
    # The clocks read at the first sample start from a line through their
    # first readings, with the drift at zero: its common part is not
    # observed, and what the start gives it stays in the scale. A later
    # clock starts at its first reading, and is NaN in the result before it.
    started = isfinite.(cc.x[1, :])
    for (i, m) in enumerate(models)
        φ, q, g = state_transition(m, τ), process_noise(m, τ), _reading_row(m)
        n = size(φ, 1)
        x0, P0 = started[i] ? _model_start(cc.x[:, i], m; zero_drift=true) : (zeros(n), zeros(n, n))
        F[i, idx(i, 2)] = 1.0
        for s in _npoly(m) .+ (1:length(_fm_taus(m)))
            F[i, idx(i, s)] = 1.0
        end
        for a in 1:n
            X̂[idx(i, a)] = x0[a]
            G[i, idx(i, a)] = g[a]
            for b in 1:n
                Φ[idx(i, a), idx(i, b)] = φ[a, b]
                Q[idx(i, a), idx(i, b)] = q[a, b]
                P[idx(i, a), idx(i, b)] = P0[a, b]
            end
        end
    end
    r = Float64[_reading_var(m) for m in models]
    X, Y, W, E, Xh = fill(NaN, N, k), fill(NaN, N, k), zeros(N, k), fill(NaN, N, k), fill(NaN, N, k)
    S = Matrix{Float64}(I, ns * k, ns * k)
    a = λ = fill(1.0 / k, k)
    for n in 1:N
        if n > 1
            X̂ = Φ * X̂
            P = Φ * P * Φ' + Q
        end
        # A clock's first reading after the start: it is the clock of most
        # weight read with it, plus a line through its first readings against
        # that clock, so that its start is in the scale's terms and the
        # reference of the readings has no part in it. The line's covariance
        # is widened 10⁵ times: those readings are still to come, and they,
        # not the line, set the clock's phase and frequency.
        for i in 1:k
            (started[i] || !isfinite(cc.x[n, i])) && continue
            d = cc.x[n:N, i:i] .- cc.x[n:N, :]
            others = [j for j in 1:k if started[j] && isfinite(d[1, j]) && count(isfinite, view(d, :, j)) > _npoly(models[i])]
            isempty(others) && continue
            j = others[argmax(a[others])]
            x0, P0 = _model_start(d[:, j], models[i]; widen=1e5)
            own = [idx(i, s) for s in eachindex(x0)]
            T = Matrix{Float64}(I, ns * k, ns * k)
            T[own, :] .= 0
            T[own[1], :] = G[j, :]
            T[own[2], :] = F[j, :]
            _npoly(models[i]) == 3 && _npoly(models[j]) == 3 && (T[own[3], idx(j, 3)] = 1.0)
            X̂ = T * X̂
            X̂[own] += x0
            P = T * P * T'
            P[own, own] += P0
            started[i] = true
        end
        use = [i for i in 1:k if started[i] && isfinite(cc.x[n, i])]
        x̂⁻ = G * X̂
        if length(use) >= 2
            # Each clock against the first one read: z = xⱼ − xₚ + vⱼ − vₚ.
            p, rest = use[1], use[2:end]
            m = length(rest)
            H = G[rest, :] .- G[p:p, :]
            R = fill(r[p], m, m)
            for (row, j) in enumerate(rest)
                R[row, row] += r[j]
            end
            z = cc.x[n, rest] .- cc.x[n, p]
            D = Symmetric(H * P * H' + R)
            K = _rdiv_cov(P * H', Matrix(D))
            X̂ = X̂ + K * (z - H * X̂)
            P = P - K * H * P
            P = (P + P') ./ 2
            # Noiseless readings leave every corrected clock the same, a
            # weighted mean of the prediction errors: a row of I − KH.
            λ = vec(mean(I - K[1:k, :] * H[:, 1:k]; dims=1))
        end
        # Reduce the phase covariance about the implicit ensemble mean of the
        # clocks read: a clock without a reading has no weight (Greenhall 2006
        # §10). Among clocks with noiseless readings the split of weight is
        # free; the gain's is reported (all the clocks, when every reading is
        # noiseless).
        if !isempty(use)
            split = length(use) >= 2 ? λ[use] : [1.0]
            aᵣ = _iem_weights(P[use, use], Bool[r[i] == 0 for i in use], split)
            a = zeros(k)
            a[use] = aᵣ !== nothing ? aᵣ : split ./ sum(split)
        end
        if algorithm === :greenhall
            S[1:k, 1:k] = I - ones(k) * a'
            P = S * P * S'
        end
        # Clock i − scale: the reading less the pivot's offset, the pivot's
        # reading less its estimated phase.
        p = pivot in use ? pivot : (isempty(use) ? 0 : use[1])
        u = p == 0 ? NaN : X̂[p] - cc.x[n, p]     # no clock of the scale is read
        ŷ = F * X̂                             # frequency: its own state plus the Markov frequency states
        for i in findall(started)
            X[n, i] = isfinite(cc.x[n, i]) ? cc.x[n, i] + u : X̂[i]
            Xh[n, i], Y[n, i] = X̂[i], ŷ[i]
        end
        n > 1 && (E[n, use] = X[n, use] .- x̂⁻[use])
        W[n, :] .= algorithm === :jones_tryon ? NaN : a
    end
    return TimescaleResult(algorithm, cc, collect((0:N-1) .* τ), X, Y, W, E, Xh)
end

# Inverse-variance weights summing to 1 over the clocks `use`: weight ∝ 1/v;
# a clock with zero variance takes all of the weight (equally among such
# clocks), and without any variance the weights are equal.
function _inverse_weights(v::Vector{Float64}, use::Vector{Int}, wmax::Float64)
    w = zeros(length(v))
    zero = [i for i in use if v[i] <= 0]
    if !isempty(zero)
        w[zero] .= 1.0
    elseif all(i -> isfinite(v[i]), use)
        w[use] .= 1.0 ./ v[use]
    else
        w[use] .= 1.0
    end
    w ./= sum(w)
    return wmax < 1 ? _cap_weights!(w, wmax) : w
end

# Start of a clock's difference against the pivot, of model `m`: the model's
# estimate at the difference's first reading from the readings of the `stop`
# samples from there (16 readings at least), which is the Kalman filter over
# them smoothed back, then carried back to sample 1 with zero added noise. The
# filter starts wide, so that it does not count its first readings twice.
function _stein_start(d::Vector{Float64}, m::ClockModel, stop::Int)
    idx = findall(isfinite, d)
    z = d[idx[1]:max(min(length(d), idx[1] + stop - 1), idx[min(length(idx), 16)])]
    x0, P0 = _model_start(z, m; widen=1e5)
    ks = kalman_smooth(kalman_filter(PhaseData(z, m.tau0), m; x0, P0))
    Φ = state_transition(m, -(idx[1] - 1) * m.tau0)
    return Φ * ks.x[1, :], Φ * ks.P[:, :, 1] * Φ'
end

# Stein's Kalman timescale (Stein 2003, eqs 4.1–7.3; Levine 2012 §XVI, KAS-2).
# Each clock is a three-state model (phase, frequency, aging; two states when no
# model has aging noise) driven by independent shocks of variance ε² = q1τ,
# η² = q2τ, α² = q3τ. Every clock i is differenced against the pivot p and the
# difference states are estimated by their own Kalman filter (5.8). The shocks
# of the pivot against the ensemble then follow from the constraint that the
# weighted shocks of the clocks sum to zero (4.1, 7.1): for each state
# ε̂ₚ = −Σᵢ aᵢ Kᵢ νᵢ with the weights a, b, c of phase, frequency and aging
# (6.2, 7.2); the other clocks are the pivot plus their differences. The weights
# are proportional to the inverse of the matching noise variance, cut at
# `wmax` as in KAS-2.
function _stein(cc::ClockComparisons, models::AbstractVector{ClockModel}, pivot::Int, wmax::Float64)
    N, k = size(cc.x)
    length(models) == k || throw(ArgumentError("timescale: $(length(models)) models for $k clocks"))
    τ = cc.tau0
    all(m -> m.tau0 ≈ τ, models) || throw(ArgumentError("timescale: the models' tau0 must be the comparisons' τ₀"))
    k * wmax >= 1 || throw(ArgumentError("timescale: wmax below 1/$k leaves weights that cannot sum to 1"))
    any(m -> _nstates(m) > _npoly(m), models) && throw(ArgumentError(
        "timescale: :stein has no flicker states; build the models with clock_model(...; flicker=false), or use :greenhall or :jones_tryon"))
    ns = any(m -> m.q3 > 0, models) ? 3 : 2
    Φ = zeros(ns, ns)                                 # element-wise: matrix literals are not trim-safe
    for d in 1:ns; Φ[d, d] = 1.0; end
    Φ[1, 2] = τ
    ns == 3 && (Φ[2, 3] = τ; Φ[1, 3] = τ^2 / 2)
    V = zeros(ns, k)                                  # shock variances per state and clock
    for (i, m) in enumerate(models)
        V[1, i], V[2, i] = (m.q1 + _fm_q1(m)) * τ, m.q2 * τ     # with the flicker FM that has no state
        ns == 3 && (V[3, i] = m.q3 * τ)
    end
    r = Float64[_reading_var(m) for m in models]
    p = pivot
    isfinite(cc.x[1, p]) || throw(ArgumentError("timescale: the pivot clock $p has no reading at the scale's first sample"))
    X, Y, W, E, Xh = fill(NaN, N, k), fill(NaN, N, k), zeros(N, k), fill(NaN, N, k), fill(NaN, N, k)
    S = zeros(ns, k)                                  # difference states, clock i − pivot
    P = [zeros(ns, ns) for _ in 1:k]
    for i in 1:k
        i == p && continue
        d = cc.x[:, i] .- cc.x[:, p]
        count(isfinite, d) > ns || continue
        md = ClockModel(τ, r[i] + r[p], (V[1, i] + V[1, p]) / τ, (V[2, i] + V[2, p]) / τ, ns == 3 ? (V[3, i] + V[3, p]) / τ : 0.0, ns == 3)
        S[:, i], P[i] = _stein_start(d, md, min(N, max(16, N ÷ 10)))
    end
    # Pivot against the scale: the scale starts as the weighted mean of the
    # clocks read at the first sample, state by state.
    members = findall(isfinite, cc.x[1, :])
    xs = zeros(ns)
    for m in 1:ns
        w = _inverse_weights(V[m, :], members, wmax)
        m == 1 && (W[1, :] = w)
        for i in members
            xs[m] -= w[i] * S[m, i]
        end
    end
    from = [something(findfirst(isfinite, view(cc.x, :, i)), N + 1) for i in 1:k]   # each clock's first reading
    for n in 1:N
        use = [i for i in 1:k if isfinite(cc.x[n, i])]
        xs⁻ = Φ * xs
        pred = Float64[xs⁻[1] + (i == p ? 0.0 : (Φ * S[:, i])[1]) for i in 1:k]
        if n > 1
            ν, K = zeros(k), zeros(ns, k)
            live = p in use
            for i in 1:k
                i == p && continue
                Pm = Φ * P[i] * Φ' + Diagonal(V[:, i] .+ V[:, p])
                if live && i in use
                    ν[i] = cc.x[n, i] - cc.x[n, p] - (Φ * S[:, i])[1]
                    K[:, i] = Pm[:, 1] ./ (Pm[1, 1] + r[i] + r[p])
                    S[:, i] = Φ * S[:, i] + K[:, i] .* ν[i]
                    P[i] = Pm - K[:, i] * Pm[1:1, :]
                    P[i] = (P[i] + P[i]') ./ 2
                else
                    S[:, i] = Φ * S[:, i]
                    P[i] = Pm
                end
            end
            xs = xs⁻
            # A clock back from a gap joins anew: its innovation spans the gap,
            # not one interval, and has no part in the pivot's shock (all of
            # them have, when every clock read is back from a gap).
            kept = [i for i in use if isfinite(cc.x[n-1, i])]
            isempty(kept) && (kept = use)
            if live && length(use) >= 2
                for m in 1:ns
                    w = _inverse_weights(V[m, :], kept, wmax)
                    xs[m] -= sum(w[i] * K[m, i] * ν[i] for i in kept)
                    m == 1 && (W[n, :] .= w)
                end
            elseif n > 1
                W[n, :] .= W[n-1, :]
            end
        end
        u = isfinite(cc.x[n, p]) ? xs[1] - cc.x[n, p] : 0.0
        for i in 1:k
            n < from[i] && continue                   # NaN before it
            Xh[n, i] = xs[1] + (i == p ? 0.0 : S[1, i])
            X[n, i] = isfinite(cc.x[n, i]) && isfinite(cc.x[n, p]) ? cc.x[n, i] + u : Xh[n, i]
            Y[n, i] = xs[2] + (i == p ? 0.0 : S[2, i])
        end
        n > 1 && (E[n, use] = X[n, use] .- pred[use])
    end
    return TimescaleResult(:stein, cc, collect((0:N-1) .* τ), X, Y, W, E, Xh)
end
