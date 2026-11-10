# Timescales from a clock ensemble. EWFA and AT1 average the clocks'
# predicted phases, with equal or learned inverse-variance weights; the
# reduced Kalman timescale runs one Kalman filter over all the clocks.

"""
    TimescaleResult

Output of [`timescale`](@ref).

$(TYPEDFIELDS)
"""
struct TimescaleResult
    "Algorithm that formed the scale (`:ewfa`, `:at1`, `:rkt`)."
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

# A clock's starting prediction-error variance from its ADEV at τ₀ over the
# first samples: the Groslambert estimate with three or more clocks, else (or
# where that comes out negative) its ADEV against the reference, which has the
# right scale but includes the reference's noise.
function _initial_variance(cc::ClockComparisons, n::Int)
    k = size(cc.x, 2)
    head = ClockComparisons(cc.x[1:n, :], cc.tau0; names=cc.names)
    σ = Float64[adev(head[i], [1]; ci=false).dev[1] for i in 1:k]
    if k >= 3
        pairs = Matrix{PhaseData{Float64}}(undef, k, k)
        for i in 1:k, j in i+1:k
            pairs[i, j] = head[i, j]
        end
        g = Float64[r.dev[1] for r in groslambert(pairs, [1])]
        σ = Float64[isfinite(g[i]) && g[i] > 0 ? g[i] : σ[i] for i in 1:k]
    end
    ok = filter(s -> isfinite(s) && s > 0, σ)
    isempty(ok) && return ones(k)
    σ = Float64[isfinite(s) && s > 0 ? s : maximum(ok) for s in σ]    # the reference itself reads zero
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

With two clocks the comparisons alone cannot tell which clock is the better
one, and AT1's weights wander toward either; cap them with `wmax` or use
EWFA.

Each clock's frequency is an exponential average of its phase steps over
`tau_y` seconds (default 100 τ₀); `tau_w` defaults to 10 `tau_y`. With
three or more clocks, AT1 starts from each clock's Groslambert ADEV over
the first 10 % of the record.

- `:rkt`, the reduced Kalman timescale in Greenhall's IEMx form
  [greenhall-2003-jones-tryon, greenhall-2006-measurement-noise](@cite):
  one Kalman filter over every clock's phase and frequency (and drift),
  with one [`ClockModel`](@ref) per clock in `models`, whose `R` is the
  clock's white PM against the reference. After each update the phase
  part of the covariance is reduced about the implicit ensemble mean, a
  weighted mean of the corrected clocks; this keeps the frequency
  estimates and minimizes each step of the scale, so the scale is quiet at
  short τ as well as long. The weights are those of that mean; `tau_y`,
  `tau_w`, `wmax` and `outlier` do not apply. The scale's frequency against
  the reference settles at a constant set by the start, which comes from a
  line through each clock's first samples
  [greenhall-2007-measurement-noise-corrections](@cite).
"""
function timescale(cc::ClockComparisons; algorithm::Symbol=:at1, tau_y::Real=100cc.tau0,
                   tau_w::Real=10tau_y, wmax::Real=1.0, outlier::Real=4.0,
                   models::Union{Nothing,AbstractVector{ClockModel}}=nothing)
    algorithm in (:ewfa, :at1, :rkt) ||
        throw(ArgumentError("timescale: algorithm must be :ewfa, :at1 or :rkt, got :$algorithm"))
    N, k = size(cc.x)
    N >= 3 || throw(ArgumentError("timescale: needs three or more samples"))
    if algorithm === :rkt
        models === nothing && throw(ArgumentError("timescale: :rkt needs models, one ClockModel per clock"))
        return _reduced_kalman(cc, models)
    end
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

function _reduced_kalman(cc::ClockComparisons, models::AbstractVector{ClockModel})
    N, k = size(cc.x)
    length(models) == k || throw(ArgumentError("timescale: $(length(models)) models for $k clocks"))
    τ = cc.tau0
    all(m -> m.tau0 ≈ τ, models) || throw(ArgumentError("timescale: the models' tau0 must be the comparisons' τ₀"))
    ns = any(m -> m.drift, models) ? 3 : 2
    # State: the k phases, then the k frequencies (then the k drifts).
    idx(i, s) = (s - 1) * k + i
    Φ = Matrix{Float64}(I, ns * k, ns * k)
    Q = zeros(ns * k, ns * k)
    X̂ = zeros(ns * k)
    P = zeros(ns * k, ns * k)
    for (i, m) in enumerate(models)
        φ, q = state_transition(m, τ), process_noise(m, τ)
        n = size(φ, 1)
        x0, P0 = _kalman_start(cc.x[:, i], τ, n)
        for a in 1:n
            X̂[idx(i, a)] = x0[a]
            for b in 1:n
                Φ[idx(i, a), idx(i, b)] = φ[a, b]
                Q[idx(i, a), idx(i, b)] = q[a, b]
                P[idx(i, a), idx(i, b)] = P0[a, b]
            end
        end
    end
    r = Float64[m.R for m in models]
    X, Y, W, E = zeros(N, k), zeros(N, k), zeros(N, k), fill(NaN, N, k)
    S = Matrix{Float64}(I, ns * k, ns * k)
    a = λ = fill(1.0 / k, k)
    for n in 1:N
        if n > 1
            X̂ = Φ * X̂
            P = Φ * P * Φ' + Q
        end
        x̂⁻ = X̂[1:k]
        use = findall(isfinite, cc.x[n, :])
        if length(use) >= 2
            # Each clock against the first one read: z = xⱼ − xₚ + vⱼ − vₚ.
            p, rest = use[1], use[2:end]
            m = length(rest)
            H = zeros(m, ns * k)
            R = fill(r[p], m, m)
            for (row, j) in enumerate(rest)
                H[row, j], H[row, p] = 1.0, -1.0
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
        S[1:k, 1:k] = I - ones(k) * a'
        P = S * P * S'
        # Clock i − scale: the reading plus reference − scale where there is
        # one, else the filter's phase.
        u = isempty(use) ? 0.0 : -sum(a[use] .* (cc.x[n, use] .- X̂[use])) / sum(a[use])
        for i in 1:k
            X[n, i] = isfinite(cc.x[n, i]) ? cc.x[n, i] + u : X̂[i]
        end
        n > 1 && (E[n, use] = X[n, use] .- x̂⁻[use])
        Y[n, :] = X̂[k+1:2k]
        W[n, :] = a
    end
    return TimescaleResult(:rkt, cc, collect((0:N-1) .* τ), X, Y, W, E)
end
