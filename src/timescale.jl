# Timescales from a clock ensemble. EWFA and AT1 average the clocks'
# predicted phases, with equal or learned inverse-variance weights; the
# reduced Kalman timescale runs one Kalman filter over all the clocks.

"""
    TimescaleResult

Output of [`timescale`](@ref).

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
    "Weight of each clock at each sample (rows sum to 1)."
    weights::Matrix{Float64}
    "Each clock's prediction error (s); `NaN` where it had no reading."
    error::Matrix{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", r::TimescaleResult)
    n, k = size(r.x)
    print(io, uppercase(String(r.algorithm)), " timescale from ", k, " clocks, ", n, " samples")
    @printf(io, "\n  final weights %s", join((@sprintf("%.3f", w) for w in r.weights[end, :]), ", "))
end

"""
    timescale_phase(r::TimescaleResult) → PhaseData

The scale against the comparisons' reference: reference − clock i plus
clock i − scale, the same for every clock read at that sample (`NaN` where
none was). Its stability is the scale's,
as far as the reference is better than it; against an ideal reference
(as from [`simulate_ensemble`](@ref) with `reference = 0`) it is exact.
"""
function timescale_phase(r::TimescaleResult)
    c = r.comparisons.x
    u = Float64[(i = findfirst(isfinite, view(c, n, :))) === nothing ? NaN : r.x[n, i] - c[n, i]
                for n in axes(c, 1)]
    return PhaseData(u, r.comparisons.tau0; source=uppercase(String(r.algorithm)))
end

"""
    implicit_mean(r::TimescaleResult, truth) → PhaseData

The implicit ensemble mean of a `:greenhall` scale, in simulation: the weighted mean of
the clocks' true phases, ``Σ aᵢ xᵢ``, with `r`'s weights and `truth` the
states from [`simulate_ensemble`](@ref) (`truth[k, 1, i]` is clock `i`'s phase
against true time). It is the scale's phase against true time without any
clock's reading noise, so its stability is the scale's own. The readout
[`timescale_phase`](@ref) at a pivot clock carries that clock's reading noise
besides; with `R = 0` they are the same up to sign (`timescale_phase` is
reference − scale) and the constant frequency offset the scale starts with. Only for `:greenhall`, whose weights define the mean.
"""
function implicit_mean(r::TimescaleResult, truth::AbstractArray{<:Real,3})
    r.algorithm === :greenhall ||
        throw(ArgumentError("implicit_mean: the implicit mean belongs to :greenhall, not :$(r.algorithm)"))
    size(truth, 1) == size(r.weights, 1) && size(truth, 3) == size(r.weights, 2) ||
        throw(ArgumentError("implicit_mean: truth is $(size(truth)), the scale has $(size(r.weights, 1)) samples and $(size(r.weights, 2)) clocks"))
    m = Float64[sum(r.weights[n, i] * truth[n, 1, i] for i in axes(r.weights, 2)) for n in axes(r.weights, 1)]
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
# from the second differences of its readings against the reference over the
# triples that are all read. The reference's own share, the mean cross product
# of two clocks' differences, is taken out; where that would leave nothing,
# the raw value stands. NaN where a clock has fewer than 8 triples.
function _time_variances(r::AbstractMatrix{Float64}, m::Int)
    N, k = size(r)
    N > 2m || return fill(NaN, k)
    D = r[1+2m:N, :] .- 2 .* r[1+m:N-m, :] .+ r[1:N-2m, :]
    raw = Float64[]
    for i in 1:k
        d = filter(isfinite, view(D, :, i))
        push!(raw, length(d) < 8 ? NaN : sum(abs2, d) / (2 * length(d)))
    end
    cross = Float64[]
    for i in 1:k, j in i+1:k
        ok = isfinite.(D[:, i]) .& isfinite.(D[:, j])
        count(ok) >= 8 && push!(cross, sum(D[ok, i] .* D[ok, j]) / (2 * count(ok)))
    end
    ref = isempty(cross) ? 0.0 : sum(cross) / length(cross)
    return Float64[v - ref > 0 ? v - ref : v for v in raw]
end

# Weiss & Weissert 1991 eq. 10: the exponential constant that best predicts a
# clock with white FM and random-walk FM whose ADEV has its minimum at τ_min.
_filter_constant(τmin::Real, τ₀::Real) = (-1 + sqrt(1 / 3 + 4 * (τmin / τ₀)^2 / 3)) / 2

# Each clock's own τ_min: the lag where its individual ADEV is smallest over
# the octaves the record supports.
function _estimate_tau_min(cc::ClockComparisons)
    N, k = size(cc.x)
    best, tmin = fill(Inf, k), fill(cc.tau0, k)
    m = 1
    while N - 2m >= 16
        v = _time_variances(cc.x, m)
        for i in 1:k
            s = sqrt(v[i]) / m
            isfinite(s) && s < best[i] && (best[i] = s; tmin[i] = m * cc.tau0)
        end
        m *= 2
    end
    return tmin
end

# Starting prediction-error variances ε_i²(τ₀) = τ₀² σ_y,i²(τ₀) (eq. 7), from
# the first max(16, N/10) readings; a clock with no usable estimate takes
# the largest of the others.
function _initial_variance(cc::ClockComparisons)
    N, k = size(cc.x)
    v = _time_variances(view(cc.x, 1:min(N, max(16, N ÷ 10)), :), 1)
    ok = filter(s -> isfinite(s) && s > 0, v)
    isempty(ok) && return ones(k)
    return Float64[isfinite(s) && s > 0 ? s : maximum(ok) for s in v]
end

# Starting frequencies against the scale: each clock's least-squares slope
# over the first max(16, N/10) readings less the mean slope, so clocks with
# frequency offsets start tracking them (the mean of the offsets is the
# scale's own).
function _initial_frequency(cc::ClockComparisons)
    N, k = size(cc.x)
    rows = 1:min(N, max(16, N ÷ 10))
    slope = zeros(k)
    for i in 1:k
        ok = isfinite.(cc.x[rows, i])
        count(ok) >= 2 || continue
        t = (rows[ok] .- 1) .* cc.tau0
        slope[i] = sum((t .- sum(t) / length(t)) .* cc.x[rows[ok], i]) / sum(abs2, t .- sum(t) / length(t))
    end
    seen = isfinite.(cc.x[1, :])
    return seen .* (slope .- sum(slope[seen]) / max(1, count(seen)))
end

# Steady-state frequency variance P of the AT2 filter (eq 16 for any τ): with
# s = σ_β²τ, the fixed point of P̂ = P + s, P = σα² P̂/(σα² + P̂), i.e.
# P² + sP − σα²s = 0. Units: σα² and s are fractional frequency squared.
_at2_steady(σα²::Float64, s::Float64) = (-s + sqrt(s^2 + 4 * σα² * s)) / 2

"""
    timescale(cc::ClockComparisons; algorithm=:at1, …) → TimescaleResult

Form a timescale from the clocks in `cc`. At each sample every clock's
phase against the scale is predicted from the last one and its estimated
frequency, and the scale is placed so that the weighted prediction errors
sum to zero. A clock without a reading keeps its prediction and drops out
of the average for that sample; a clock joins at its first reading. Each
clock's frequency starts at its least-squares slope over the first 10 % of
the record, less the mean slope.

- `:ewfa`, equal-weights frequency averaging, section 3.1 of
  [kettering-2023-clock-ensemble](@cite): every clock weighs the same, and the
  scale's frequency is the mean of the clocks' frequency steps over each
  interval, so with every clock read at every sample the scale is the plain
  mean of the clocks' phases. Gaps are not in the thesis: here a clock without
  a reading keeps its prediction (its last frequency step carried on) and
  drops out of the average for that sample, so the scale does not jump when a
  clock drops out, as the plain mean would.
- `:at1`, the NIST scale as Weiss & Weissert describe it
  [weiss-weissert-1991-at2](@cite) (equation numbers are the paper's):
  weights ``w_i = ε_x²/ε_i²`` (3), with ``ε_x² = (Σ 1/ε_i²)⁻¹`` (4), where
  ``ε_i²`` is an exponential average (7) of the squared prediction error
  estimate ``ε̂_i = |x̂_i − x_i| + K_i`` (5), the bias
  ``K_i = 2ε_x²/√(2πε_i²)`` (6) making up for each clock being part of the
  scale. Each clock's frequency is filtered (9) with its own constant from
  its τ_min (10).
- `:at2`, AT1 plus a frequency variance
  [weiss-weissert-1991-at2](@cite): the weights are AT1's, and each clock's
  frequency ``y`` is a Kalman filter on its average frequency ``ŷ`` (8): the
  variance ``P`` grows by the random-walk FM impulse ``σ_β²τ`` (13), then
  ``y ← (σ_α² y + P̂ ŷ)/(σ_α² + P̂)`` (14) and ``P ← σ_α² P̂/(σ_α² + P̂)`` (15), with the
  white FM ``σ_α² = ε_i²/τ²`` from the prediction error, so the filter
  adapts as ``ε_i`` does. ``σ_β²τ`` is fixed per clock from its τ_min by the
  steady state ``m(m+1) = σ_α²/(σ_β²τ)`` (18), the starting ``σ_α²`` being the
  starting variance, and ``P`` starts at its steady state (16). In steady state
  this is AT1's filter (17). With `steps=true` the paper's frequency-step
  detection (19–22) runs too, see below. Takes the AT1 keywords.
- `:greenhall`, Greenhall's Kalman timescale
  [greenhall-2003-jones-tryon, greenhall-2006-measurement-noise](@cite):
  one Kalman filter over every clock's phase and frequency (and drift),
  with one [`ClockModel`](@ref) per clock in `models`, whose `R` is the
  clock's white PM against the reference. After each update the phase
  part of the covariance is reduced about the implicit ensemble mean, a
  weighted mean of the corrected clocks; this keeps the frequency
  estimates and minimizes each step of the scale, so the scale is quiet at
  short τ as well as long. The weights are those of that mean. The scale's
  frequency against the reference settles at a constant set by the start,
  which comes from a line through each clock's first samples
  [greenhall-2007-measurement-noise-corrections](@cite). With every `R = 0`
  this is the reduced Kalman timescale of the 2003 paper; with `R > 0` it is
  IEMx (2006), the same filter. Readings with noise leave the corrected clocks
  unequal and the implicit ensemble mean unobservable, so the scale is
  realized at clock `pivot` (default 1; the first clock read when it has no
  reading): that clock's reading less its estimated offset from the implicit
  mean. It carries that clock's own reading noise as white PM, so it is not
  quieter than a weighted readout (three equal clocks: √3 noisier at short
  τ), and pivots differ by that noise; at `R = 0` every pivot gives the same
  scale. While the pivot has no reading the first clock read stands in; a
  pivot that starts late steps the scale at its first reading. Models with
  flicker FM states bring them into the filter: at `R = 0` that is the
  algorithm of [davis-greenhall-stacey-2005-flicker](@cite) (its eqs 3.1–3.6
  and 4.1), and with `R > 0` the 2006 paper's reduction covers the same
  state (its §2: states that reach the readings only through the phases).
  The frequency reported for a clock then includes them. Flicker PM states
  are measurement noise with memory: they enter the readings, not the phases
  the scale is formed from, so they leave the weights alone. That is outside
  the 2006 paper, whose measurement noise is white; the reduction only needs
  the readings to be differences of clocks, which they still are.

- `:jones_tryon`, the natural Jones–Tryon timescale
  [jones-tryon-1983-estimating-time, greenhall-2003-jones-tryon](@cite): the same
  filter as `:greenhall` without the covariance reduction. The scale is a
  clock's reading less its estimated phase (JT83 §4, after eq. 4.17;
  Greenhall 2003 eq. 10); with noiseless readings it is the same at every
  clock. The frequency estimates are the same as `:greenhall`'s; the scale is
  slightly noisier at short τ (the reduction is what minimizes its steps:
  about 4% at τ = 1 on four white-FM clocks, ten seeds, with single records
  up to 40% higher). Its weights are `NaN`: no weighted mean defines it. `R`
  is the white PM of a clock's reading against the reference, so every
  difference carries the reference's noise (correlated), where JT83 takes
  independent errors in each difference. JT83's error
  detection, step corrections and drift state are not part of it.
- `:stein`, Stein's Kalman timescale [stein-2003-demystified](@cite), the
  KAS-2 algorithm [levine-2012-statistical-modeling, erickson-2011-symmetricom-time-scale](@cite),
  with one [`ClockModel`](@ref) per clock in `models`. Each clock has phase,
  frequency and aging states (no aging when no model has `q3 > 0`) driven by
  independent shocks of variance `q1 τ`, `q2 τ` and `q3 τ`; `R` is the clock's
  reading noise. The clocks are differenced against clock `pivot` and each
  difference gets its own Kalman filter (Stein 2003 eq. 5.8). The pivot's own
  shocks against the ensemble follow from the constraints that the weighted shocks
  of all the clocks sum to zero (4.1, 7.1), one set of weights for each state, and
  the other clocks are the pivot plus their difference estimates (6.3, 7.3). The
  weights are inversely proportional to the variance of the matching noise
  (Levine §XVI; a clock with no such noise takes all of it); `wmax` caps them as
  in KAS-2. The weights returned are the phase weights. The scale is realized at
  the pivot as for `:greenhall`. Departures from the paper: the shock variances
  are `qτ` (no cross terms between the states, as in the paper's shock model,
  where `:greenhall` integrates them), a clock without a reading in a step is
  predicted from its difference state and the weights run over the clocks read,
  and the paper's independent difference filters ignore the correlation of
  their innovations through the pivot, so the scale differs a little with the
  pivot. With white FM clocks it matches `:greenhall` to a few percent at every τ.

AT1 keywords:

- `tau_min`: each clock's τ_min in seconds, a number or one per clock; the
  paper takes it from the clock's characterization. Default: the lag of the
  smallest ADEV of each clock, estimated from the record (the reference's
  share taken out).
- `memory`: the time constant of the variance average in seconds, the
  paper's ``N_τ`` (20 days for cesium clocks, and the default; at most a
  quarter of the record).
- `bias`: `:paper` (default) for ``K_i`` above, or `:exact` for the squared
  error divided by ``1 − w_i``, which is exact when the clocks' own errors
  are independent and the weights optimal.
- `steps` (`:at2`, default `false`): detect a clock's frequency step and rerun
  the scale without it, in one of two modes. `:best` (`true`) replaces the
  scale from the step on, the best estimate of a recorded history. `:real`
  keeps every value already issued, as the paper does (it cannot change
  history and steers forward): the rerun only corrects what follows, continuing from the
  last issued phase, so a step costs the samples before it was caught.
  After each sample, for windows of 2 to `lmax − 1`
  intervals (`lmax` defaults to the clock's τ_min in samples, at least 4),
  the clock's mean frequency against the scale over the window is compared
  with its frequency estimate before the window (eqs 19–20). A difference
  above `step_threshold` (default 5) standard deviations is a step at the
  window's start: the scale is rerun from there with that clock's weight
  zero until τ_min later and its frequency variance raised, and it comes back
  with ``ε_i²`` multiplied by `step_inflate` (default 4, σ_α doubled), the
  paper's "increased by a fixed amount". Departures: eq 20 as printed mixes
  units (τ_min/L is a time), so the deviation is the variance of the
  difference, in fractional frequency squared,
  ``P_{−L} + (σ_{α,i}² + σ_{αx}²)/s + (σ_{β,i}²τ + σ_{βx}²τ)(s+1)(2s+1)/(6s)``
  for a window of ``s`` intervals: the estimate's variance, the white FM of
  the mean over the window (clock and scale, eq 21), and the frequency random
  walk over it (the mean of a random walk of impulse σ_β²τ per interval over
  ``s`` intervals has variance σ_β²τ (s+1)(2s+1)/(6s), about s/3; scale term
  as eq 22). The paper's ``(P_{-L}+P_{-2})/2`` is ``P_{-L}`` alone here, and
  windows start at two intervals, my reading of its two-back start (one
  interval would flag more small phase jumps; a large lone jump is still
  flagged). The clock-minus-scale variance is
  less than the sum of the two when the clock is in the scale, so the sum is
  conservative, as in the paper. A clock is
  not tested over windows that reach into its exclusion. The paper's constant is 4, which gives about one false step per
  3×10⁴ clock-samples on the noisier clocks of a step-free ensemble (53 in 24
  runs of four clocks and 16384 samples, none at 5 or 6); a false step
  excludes the clock for its τ_min and raises its ε², which shifts the
  weights for good. At 5, a frequency step of 5×10⁻¹¹ (several times the
  white-FM noise per interval) is caught after about 6 samples and 2×10⁻¹¹
  after about 50; one of 10⁻¹¹ is not. Pass `step_threshold=4` for the paper's value.
- `wmax` (default 1, no cap): the largest weight any clock may have.
- `outlier` (default `Inf`, off): a clock whose prediction error passes
  `outlier` ε_i is left out of that sample and of its variance; if that
  would leave none, all stay.

The weights are optimal for the interval of the readings: ε_i² is the
clock's prediction error over one interval, so for τ well beyond it, and
with clocks of different long-term behaviour, the scale can be worse than
the best clock. With two clocks the comparisons alone cannot tell which is
better, and the weights wander toward either; cap them with `wmax` or use
EWFA. The frequency drift term of the paper's prediction (1) is not used.
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
    N, k = size(cc.x)
    N >= 3 || throw(ArgumentError("timescale: needs three or more samples"))
    if algorithm === :stein
        models === nothing && throw(ArgumentError("timescale: :stein needs models, one ClockModel per clock"))
        1 <= pivot <= k || throw(ArgumentError("timescale: pivot $pivot is not one of the $k clocks"))
        return _stein(cc, models, Int(pivot), Float64(wmax))
    end
    if algorithm === :greenhall || algorithm === :jones_tryon
        models === nothing && throw(ArgumentError("timescale: :$algorithm needs models, one ClockModel per clock"))
        1 <= pivot <= k || throw(ArgumentError("timescale: pivot $pivot is not one of the $k clocks"))
        return _reduced_kalman(cc, models, Int(pivot), algorithm)
    end
    mode = steps === true ? :best : steps === false ? :none : steps
    mode in (:none, :best, :real) ||
        throw(ArgumentError("timescale: steps must be false, :best or :real, got :$steps"))
    mode !== :none && algorithm !== :at2 && throw(ArgumentError("timescale: steps needs algorithm :at2"))
    k * wmax >= 1 || throw(ArgumentError("timescale: wmax below 1/$k leaves weights that cannot sum to 1"))
    at2 = algorithm === :at2
    at1 = algorithm === :at1 || at2
    τ = cc.tau0
    τmin = tau_min === nothing ? (at1 ? _estimate_tau_min(cc) : fill(τ, k)) :
           tau_min isa Real ? fill(Float64(tau_min), k) : Float64.(tau_min)
    length(τmin) == k || throw(ArgumentError("timescale: $(length(τmin)) tau_min values for $k clocks"))
    mf = at1 ? Float64[_filter_constant(t, τ) for t in τmin] : zeros(k)
    Nτ = max(1.0, min(memory, N * τ / 4) / τ)
    r = cc.x
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    ε² = at1 ? _initial_variance(cc) : ones(k)
    # Start: the scale is the plain mean of the first readings.
    ok = isfinite.(r[1, :])
    any(ok) || throw(ArgumentError("timescale: no clock has a first reading"))
    u = -sum(r[1, ok]) / count(ok)
    X[1, :] = r[1, :] .+ u
    seen = copy(ok)             # a clock joins the scale at its first reading
    y = _initial_frequency(cc)
    # AT2: the random-walk FM impulse σ_β²τ of each clock from its τ_min through
    # the steady state m(m+1) = σα²/(σ_β²τ) (eq 18), and the frequency variance P
    # at its steady state (eq 16), both from the starting σα² = ε²/τ².
    σβτ = ε² ./ τ^2 ./ (mf .* (mf .+ 1))
    Pf = Float64[_at2_steady(ε²[i] / τ^2, σβτ[i]) for i in 1:k]
    W[1, :] = ok ./ count(ok)
    # Step detector state (AT2, eqs 19-22): the state after each sample, kept
    # to rerun the scale from a detected step; `blocked[i]` is the last sample
    # clock i is kept out of the scale; no detection until `replay`.
    Lmax = [lmax === nothing ? max(4, round(Int, t / τ)) : Int(lmax) for t in τmin]
    ε²h, yh, Pfh = zeros(N, k), zeros(N, k), zeros(N, k)
    blocked, replay = zeros(Int, k), 0
    # :real keeps what was issued: rows are copied out once, after a rerun
    # continued from the last issued phase (the correction steers forward).
    real = mode === :real
    Xo, Yo, Wo, Eo = real ? (zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)) : (X, Y, W, E)
    issued, pending, shift = 0, false, 0.0
    n = 2
    while n <= N
        if real && n - 1 > issued
            if pending     # phase continuity with the last row issued
                j = findall(isfinite, view(r, issued, :))
                shift = isempty(j) ? 0.0 : sum(Xo[issued, j] .- X[issued, j]) / length(j)
                pending = false
            end
            Xo[n-1, :], Yo[n-1, :], Wo[n-1, :], Eo[n-1, :] = X[n-1, :] .+ shift, Y[n-1, :], W[n-1, :], E[n-1, :]
            issued = n - 1
        end
        ε²h[n-1, :], yh[n-1, :], Pfh[n-1, :] = ε², y, Pf
        for i in 1:k
            n == blocked[i] + 1 && (ε²[i] *= step_inflate)   # back in, with σα raised
        end
        x̂ = X[n-1, :] .+ y .* τ
        use = isfinite.(r[n, :]) .& seen .& (n .> blocked)
        if at1 && isfinite(outlier) && n > 2 && any(use)
            # The median prediction error places the scale for the test, so
            # that one clock's jump cannot carry it.
            d = x̂ .- r[n, :]
            keep = use .& (abs.(d .- _middle(d[use])) .<= outlier .* sqrt.(ε²))
            any(keep) && (use = keep)
        end
        if !any(use)
            X[n, :] = x̂
            n += 1
            continue
        end
        wu = at1 ? use ./ ε² : Float64.(use)
        wu ./= sum(wu)
        εx² = at1 ? 1 / sum(1 ./ ε²[use]) : 0.0
        _cap_weights!(wu, max(Float64(wmax), 1 / count(use)))
        u = sum(wu[use] .* (x̂[use] .- r[n, use]))
        for i in 1:k
            if isfinite(r[n, i]) && !seen[i]
                X[n, i] = r[n, i] + u
                seen[i] = true
            elseif isfinite(r[n, i])
                X[n, i] = r[n, i] + u
                E[n, i] = X[n, i] - x̂[i]
                if at1 && use[i]
                    ε̂² = bias === :paper ? (abs(E[n, i]) + 2εx² / sqrt(2π * ε²[i]))^2 :
                         E[n, i]^2 / (1 - min(wu[i], 0.9))
                    ε²[i] = (ε̂² + Nτ * ε²[i]) / (1 + Nτ)
                end
                # a clock the outlier test left out keeps its frequency
                if use[i]
                    ŷ = (X[n, i] - X[n-1, i]) / τ
                    if at2          # Kalman update of y (eqs 13-15), σα² = ε²/τ²
                        σα² = ε²[i] / τ^2
                        P̂ = Pf[i] + σβτ[i]
                        y[i] = (σα² * y[i] + P̂ * ŷ) / (σα² + P̂)
                        Pf[i] = σα² * P̂ / (σα² + P̂)
                    else
                        y[i] = (mf[i] * y[i] + ŷ) / (mf[i] + 1)
                    end
                end
            else
                X[n, i] = x̂[i]
            end
        end
        Y[n, :] = y
        W[n, :] = wu
        if mode !== :none && n > replay
            n0, i0 = _at2_step(X, yh, Pfh, ε², σβτ, use, blocked, Lmax, n, τ, step_threshold)
            if i0 > 0       # a frequency step: rerun from its start without that clock
                blocked[i0] = n0 + max(1, round(Int, τmin[i0] / τ))
                replay = n
                pending = real
                ε², y, Pf = ε²h[n0, :], yh[n0, :], Pfh[n0, :]
                Pf[i0] = ε²[i0] / τ^2          # P raised at the step
                for i in 1:k
                    seen[i] = any(isfinite, view(r, 1:n0, i))
                end
                n = n0 + 1
                continue
            end
        end
        n += 1
    end
    if real && N > issued       # the last row, when no rerun is pending for it
        if pending
            j = findall(isfinite, view(r, issued, :))
            shift = isempty(j) ? 0.0 : sum(Xo[issued, j] .- X[issued, j]) / length(j)
        end
        for m in issued+1:N
            Xo[m, :], Yo[m, :], Wo[m, :], Eo[m, :] = X[m, :] .+ shift, Y[m, :], W[m, :], E[m, :]
        end
    end
    return TimescaleResult(algorithm, cc, collect((0:N-1) .* τ), Xo, Yo, Wo, Eo)
end

# AT2 frequency-step test (Weiss & Weissert 1991, eqs 19-22), run on the clocks
# kept at sample n: for each window of s = 2 … Lmax−1 intervals, the clock's mean
# frequency against the scale over the window, y_avg = (x[n] − x[n−s])/(sτ),
# is compared with the frequency estimate y₋ₗ before it. Eq 20 as printed
# mixes units; this is the variance of the difference, in fractional frequency²:
#   P₋ₗ + (σα,i² + σαx²)/s + (σβτ,i + σβτx)(s+1)(2s+1)/(6s)
# the estimate's variance, the white-FM noise of the mean over s intervals (the
# clock's and the scale's), and the random walk of the frequency over the window
# (the mean of a random walk of impulse σβτ per interval over s intervals has
# variance σβτ (s+1)(2s+1)/(6s)). σαx² and σβτx combine the clocks' levels as in eqs 21-22,
# with σβτ the paper's σβ²τ. Returns (start sample, clock) of the window with
# the largest ratio above `thr`, or (0, 0).
function _at2_step(X, yh, Pfh, ε², σβτ, use, blocked, Lmax, n, τ, thr)
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
        for s in 2:min(Lmax[i] - 1, n - max(1, blocked[i] + 1))
            n0 = n - s
            var = Pfh[n0, i] + (σα² + 1 / σαx⁻) / s + (σβτ[i] + 1 / σβx⁻) * (s + 1) * (2s + 1) / (6s)
            z = abs((X[n, i] - X[n0, i]) / (s * τ) - yh[n0, i]) / sqrt(var)
            z > best && ((best, n0b, ib) = (z, n0, i))
        end
    end
    return n0b, ib
end

# IEM weights a (Σa = 1) minimizing aᵀ Pxx a: Pxx⁻¹1 / 1ᵀPxx⁻¹1; `nothing`
# when Pxx is singular, as noiseless readings make it.
function _iem_weights(Pxx::Matrix{Float64})
    s = maximum(diag(Pxx))
    s > 0 || return nothing
    C = cholesky(Symmetric(Pxx ./ s); check=false)
    (issuccess(C) && minimum(diag(C.U)) > 1e-6) || return nothing
    a = C \ ones(size(Pxx, 1))
    return a ./ sum(a)
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
    G = zeros(k, ns * k)                      # each clock's reading in state terms
    for (i, m) in enumerate(models)
        φ, q, g = state_transition(m, τ), process_noise(m, τ), _reading_row(m)
        n = size(φ, 1)
        x0, P0 = _model_start(cc.x[:, i], m)
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
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    S = Matrix{Float64}(I, ns * k, ns * k)
    a = λ = fill(1.0 / k, k)
    for n in 1:N
        if n > 1
            X̂ = Φ * X̂
            P = Φ * P * Φ' + Q
        end
        x̂⁻ = G * X̂
        use = findall(isfinite, cc.x[n, :])
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
        # Reduce the phase covariance about the implicit ensemble mean; any
        # weights do when it is singular, and the gain's are reported.
        aᵣ = _iem_weights(P[1:k, 1:k])
        a = aᵣ !== nothing ? aᵣ : length(use) >= 2 ? λ : a
        if algorithm === :greenhall
            S[1:k, 1:k] = I - ones(k) * a'
            P = S * P * S'
        end
        # Clock i − scale: the reading less the pivot's offset, the pivot's
        # reading less its estimated phase.
        p = pivot in use ? pivot : (isempty(use) ? 0 : use[1])
        u = p == 0 ? 0.0 : X̂[p] - cc.x[n, p]
        for i in 1:k
            X[n, i] = isfinite(cc.x[n, i]) ? cc.x[n, i] + u : X̂[i]
        end
        n > 1 && (E[n, use] = X[n, use] .- x̂⁻[use])
        for (i, m) in enumerate(models)       # frequency: its own state plus the Markov frequency states
            Y[n, i] = X̂[k+i]
            for s in _npoly(m) .+ (1:length(_fm_taus(m)))
                Y[n, i] += X̂[idx(i, s)]
            end
        end
        W[n, :] .= algorithm === :jones_tryon ? NaN : a
    end
    return TimescaleResult(algorithm, cc, collect((0:N-1) .* τ), X, Y, W, E)
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
    any(m -> m.flicker_fm > 0 || m.flicker_pm > 0, models) && throw(ArgumentError(
        "timescale: :stein has no flicker states; pass models without flicker, or use :greenhall or :jones_tryon"))
    ns = any(m -> m.q3 > 0, models) ? 3 : 2
    Φ = zeros(ns, ns)                                 # element-wise: matrix literals are not trim-safe
    for d in 1:ns; Φ[d, d] = 1.0; end
    Φ[1, 2] = τ
    ns == 3 && (Φ[2, 3] = τ; Φ[1, 3] = τ^2 / 2)
    V = zeros(ns, k)                                  # shock variances per state and clock
    for (i, m) in enumerate(models)
        V[1, i], V[2, i] = m.q1 * τ, m.q2 * τ
        ns == 3 && (V[3, i] = m.q3 * τ)
    end
    r = Float64[m.R for m in models]
    p = pivot
    isfinite(cc.x[1, p]) || throw(ArgumentError("timescale: the pivot clock $p has no first reading"))
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    # Pivot against the scale: the start scale is the plain mean of the first readings.
    ok = isfinite.(cc.x[1, :])
    xs = zeros(ns)
    xs[1], xs[2] = cc.x[1, p] - sum(cc.x[1, ok]) / count(ok), _initial_frequency(cc)[p]
    S = zeros(ns, k)                                  # difference states, clock i − pivot
    P = [zeros(ns, ns) for _ in 1:k]
    for i in 1:k
        i == p && continue
        d = cc.x[:, i] .- cc.x[:, p]
        count(isfinite, d) > ns || continue
        S[:, i], P[i] = _kalman_start(d, τ, ns)
    end
    W .= 1.0 / k
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
            if live && length(use) >= 2
                for m in 1:ns
                    w = _inverse_weights(V[m, :], use, wmax)
                    xs[m] -= sum(w[i] * K[m, i] * ν[i] for i in use)
                    m == 1 && (W[n, :] .= w)
                end
            elseif n > 1
                W[n, :] .= W[n-1, :]
            end
        end
        u = isfinite(cc.x[n, p]) ? xs[1] - cc.x[n, p] : 0.0
        for i in 1:k
            xi = xs[1] + (i == p ? 0.0 : S[1, i])
            X[n, i] = isfinite(cc.x[n, i]) && isfinite(cc.x[n, p]) ? cc.x[n, i] + u : xi
            Y[n, i] = xs[2] + (i == p ? 0.0 : S[2, i])
        end
        n > 1 && (E[n, use] = X[n, use] .- pred[use])
    end
    return TimescaleResult(:stein, cc, collect((0:N-1) .* τ), X, Y, W, E)
end
