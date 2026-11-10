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

"""
    timescale(cc::ClockComparisons; algorithm=:at1, …) → TimescaleResult

Form a timescale from the clocks in `cc`. At each sample every clock's
phase against the scale is predicted from the last one and its estimated
frequency, and the scale is placed so that the weighted prediction errors
sum to zero. A clock without a reading keeps its prediction and drops out
of the average for that sample; a clock joins at its first reading. Each
clock's frequency starts at its least-squares slope over the first 10 % of
the record, less the mean slope.

- `:ewfa`, equal-weight frequency averaging: every clock weighs the same.
  Each clock's frequency is an exponential average of its phase steps
  over `tau_y` seconds (default 100 τ₀). With every clock read at every
  sample the scale is the plain mean of the clocks' phases.
- `:at1`, the NIST scale as Weiss & Weissert describe it
  [weiss-weissert-1991-at2](@cite) (equation numbers are the paper's):
  weights ``w_i = ε_x²/ε_i²`` (3), with ``ε_x² = (Σ 1/ε_i²)⁻¹`` (4), where
  ``ε_i²`` is an exponential average (7) of the squared prediction error
  estimate ``ε̂_i = |x̂_i − x_i| + K_i`` (5), the bias
  ``K_i = 2ε_x²/√(2πε_i²)`` (6) making up for each clock being part of the
  scale. Each clock's frequency is filtered (9) with its own constant from
  its τ_min (10).
- `:rkt`, the reduced Kalman timescale in Greenhall's IEMx form
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
  [greenhall-2007-measurement-noise-corrections](@cite).

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
function timescale(cc::ClockComparisons; algorithm::Symbol=:at1, tau_y::Real=100cc.tau0,
                   tau_min::Union{Nothing,Real,AbstractVector{<:Real}}=nothing,
                   memory::Real=20 * 86400.0, bias::Symbol=:paper, wmax::Real=1.0,
                   outlier::Real=Inf, models::Union{Nothing,AbstractVector{ClockModel}}=nothing)
    algorithm in (:ewfa, :at1, :rkt) ||
        throw(ArgumentError("timescale: algorithm must be :ewfa, :at1 or :rkt, got :$algorithm"))
    bias in (:paper, :exact) || throw(ArgumentError("timescale: bias must be :paper or :exact, got :$bias"))
    N, k = size(cc.x)
    N >= 3 || throw(ArgumentError("timescale: needs three or more samples"))
    if algorithm === :rkt
        models === nothing && throw(ArgumentError("timescale: :rkt needs models, one ClockModel per clock"))
        return _reduced_kalman(cc, models)
    end
    k * wmax >= 1 || throw(ArgumentError("timescale: wmax below 1/$k leaves weights that cannot sum to 1"))
    at1 = algorithm === :at1
    τ = cc.tau0
    τmin = tau_min === nothing ? (at1 ? _estimate_tau_min(cc) : fill(τ, k)) :
           tau_min isa Real ? fill(Float64(tau_min), k) : Float64.(tau_min)
    length(τmin) == k || throw(ArgumentError("timescale: $(length(τmin)) tau_min values for $k clocks"))
    mf = at1 ? Float64[_filter_constant(t, τ) for t in τmin] : fill(tau_y / τ, k)
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
    W[1, :] = ok ./ count(ok)
    for n in 2:N
        x̂ = X[n-1, :] .+ y .* τ
        use = isfinite.(r[n, :]) .& seen
        if at1 && isfinite(outlier) && n > 2 && any(use)
            # The median prediction error places the scale for the test, so
            # that one clock's jump cannot carry it.
            d = x̂ .- r[n, :]
            keep = use .& (abs.(d .- median(d[use])) .<= outlier .* sqrt.(ε²))
            any(keep) && (use = keep)
        end
        if !any(use)
            X[n, :] = x̂
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
                y[i] = (mf[i] * y[i] + (X[n, i] - X[n-1, i]) / τ) / (mf[i] + 1)
                if at1 && use[i]
                    ε̂² = bias === :paper ? (abs(E[n, i]) + 2εx² / sqrt(2π * ε²[i]))^2 :
                         E[n, i]^2 / (1 - min(wu[i], 0.9))
                    ε²[i] = (ε̂² + Nτ * ε²[i]) / (1 + Nτ)
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
