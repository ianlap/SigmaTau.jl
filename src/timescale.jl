# Timescales from a clock ensemble: a weighted average of the clocks'
# predicted phases. EWFA weighs every clock the same; AT1 weighs by inverse
# prediction-error variance, learned as it runs.

"""
    TimescaleResult

Output of [`timescale`](@ref).

$(TYPEDFIELDS)
"""
struct TimescaleResult
    "Algorithm that formed the scale (`:ewfa`, `:at1`)."
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
clock i − scale, the same for every clock. Its stability is the scale's,
as far as the reference is better than it; against an ideal reference
(as from [`simulate_ensemble`](@ref) with `reference = 0`) it is exact.
"""
timescale_phase(r::TimescaleResult) = PhaseData(r.x[:, 1] .- r.comparisons.x[:, 1], r.comparisons.tau0;
                                                source=uppercase(String(r.algorithm)))

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

# A clock's starting prediction-error variance from the Groslambert ADEV of
# the first samples at τ₀ (three or more clocks), else the same for all.
function _initial_variance(cc::ClockComparisons, n::Int)
    k = size(cc.x, 2)
    k >= 3 || return ones(k)
    head = ClockComparisons(cc.x[1:n, :], cc.tau0; names=cc.names)
    pairs = Matrix{PhaseData{Float64}}(undef, k, k)
    for i in 1:k, j in i+1:k
        pairs[i, j] = head[i, j]
    end
    σ = [r.dev[1] for r in groslambert(pairs, [1])]
    all(s -> isfinite(s) && s > 0, σ) || return ones(k)
    return (σ .* cc.tau0) .^ 2
end

"""
    timescale(cc::ClockComparisons; algorithm=:at1, tau_y, tau_w, wmax=1, outlier=4) → TimescaleResult

Form a timescale from the clocks in `cc`. At each sample every clock's
phase against the scale is predicted from the last one and its estimated
frequency, and the scale is placed so that the weighted prediction errors
sum to zero. A clock without a reading keeps its prediction and drops out
of the average for that sample.

- `:ewfa`, equal-weight frequency averaging: every clock weighs the same.
- `:at1`, the NIST real-time scale (Allan, Gray & Machlan 1974; Weiss,
  Allan & Peppler 1989): weights ∝ 1/σᵢ², with σᵢ² an exponential average
  of the clock's squared prediction error over `tau_w` seconds, corrected
  by 1/(1 − wᵢ) for the clock's own pull on the scale. Weights are capped
  at `wmax`, and a clock whose error passes `outlier` σᵢ is left out of
  that sample and of its variance.

Each clock's frequency is an exponential average of its phase steps over
`tau_y` seconds (default 100 τ₀); `tau_w` defaults to 10 `tau_y`. With
three or more clocks, AT1 starts from each clock's Groslambert ADEV over
the first 10 % of the record.
"""
function timescale(cc::ClockComparisons; algorithm::Symbol=:at1, tau_y::Real=100cc.tau0,
                   tau_w::Real=10tau_y, wmax::Real=1.0, outlier::Real=4.0)
    algorithm in (:ewfa, :at1) || throw(ArgumentError("timescale: algorithm must be :ewfa or :at1, got :$algorithm"))
    N, k = size(cc.x)
    N >= 3 || throw(ArgumentError("timescale: needs three or more samples"))
    k * wmax >= 1 || throw(ArgumentError("timescale: wmax below 1/$k leaves weights that cannot sum to 1"))
    τ = cc.tau0
    my, mw = tau_y / τ, tau_w / τ
    r = cc.x
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    σ² = algorithm === :at1 ? _initial_variance(cc, max(16, N ÷ 10)) : ones(k)
    w = fill(1.0 / k, k)
    # Start: the scale is the plain mean of the first readings.
    ok = isfinite.(r[1, :])
    any(ok) || throw(ArgumentError("timescale: no clock has a first reading"))
    u = -sum(r[1, ok]) / count(ok)
    X[1, :] = r[1, :] .+ u
    y = zeros(k)
    W[1, :] = w
    for n in 2:N
        algorithm === :at1 && (w = 1 ./ σ²)
        x̂ = X[n-1, :] .+ y .* τ
        use = isfinite.(r[n, :])
        if algorithm === :at1 && n > 2
            # A first pass on every reading sets the scale for the outlier test.
            u₀ = sum(w[use] .* (x̂[use] .- r[n, use])) / sum(w[use])   # uncapped is close enough here
            for i in findall(use)
                abs(r[n, i] + u₀ - x̂[i]) > outlier * sqrt(σ²[i]) && (use[i] = false)
            end
        end
        if !any(use)
            X[n, :] = x̂
            W[n, :] .= 0
            continue
        end
        wu = w .* use
        wu ./= sum(wu)
        _cap_weights!(wu, max(Float64(wmax), 1 / count(use)))
        u = sum(wu[use] .* (x̂[use] .- r[n, use]))
        for i in 1:k
            if isfinite(r[n, i])
                X[n, i] = r[n, i] + u
                E[n, i] = X[n, i] - x̂[i]
                y[i] = (my * y[i] + (X[n, i] - X[n-1, i]) / τ) / (my + 1)
                if algorithm === :at1 && use[i]
                    σ²[i] = (mw * σ²[i] + E[n, i]^2 / (1 - min(wu[i], 0.9))) / (mw + 1)
                end
            else
                X[n, i] = x̂[i]
            end
        end
        Y[n, :] = y
        W[n, :] = wu
    end
    return TimescaleResult(algorithm, cc, collect((0:N-1) .* τ), X, Y, W, E)
end
