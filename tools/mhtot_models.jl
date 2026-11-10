# tools/mhtot_models.jl — exact bias and EDF of MHTOTDEV: the source of
# `_MHTOT_B`, `_MHTOT_B_SMALL`, `_coeff_mhtot` and `_MHTOT_WHPM_EDF_SMALL` in
# src/edf.jl. No random numbers.
#
#   julia -t auto tools/mhtot_models.jl        # about five minutes
#
# MHTOTVAR takes each window of 4m phase samples, removes its least-squares
# quadratic, mirrors it at both ends (12m samples) and averages the squared
# m-averaged third difference over the 8m + 1 positions. That is a quadratic
# form w'Qw which annihilates quadratics, so E[MHTOTVAR] = tr(QΣ) depends on m
# and the noise type only. Over a record the statistic is x'Q̄x with Q̄ the
# mean of the window forms, and for Gaussian noise
# EDF = tr(Q̄Σ)² / tr((Q̄Σ)²). The matrices are built from the definition,
# independently of the package kernel; test/stab/mhtotdev_models.jl holds the
# kernel to the same numbers.
using LinearAlgebra, Printf

const NOISES = ((2, "white PM"), (1, "flicker PM"), (0, "white FM"), (-1, "flicker FM"),
                (-2, "random-walk FM"), (-3, "flicker-walk FM"), (-4, "random-run FM"))

# Third difference of m-sample box sums at every start of a `len`-sample series.
function boxes(len, m)
    D = zeros(len - 4m + 1, len)
    for j in axes(D, 1)
        D[j, j:j+m-1] .= 1; D[j, j+m:j+2m-1] .= -3; D[j, j+2m:j+3m-1] .= 3; D[j, j+3m:j+4m-1] .= -1
    end
    return D
end

# M with M'M the window's form: detrend, mirror, difference, normalize.
function window_op(m)
    L = 4m
    t = range(-1, 1, length=L)
    P = I - (q = Matrix(qr(hcat(ones(L), t, t .^ 2)).Q)[:, 1:3]; q * q')
    E = zeros(3L, L)
    for k in 1:L; E[k, L+1-k] = 1; E[L+k, k] = 1; E[2L+k, L+1-k] = 1; end
    return boxes(3L, m) * (E * P) ./ sqrt((8m + 1) * 6.0 * m^4)
end

# F with Σ = F·F' for n consecutive phase samples: white noise, or flicker by
# Kasdin's half-order filter over a 4n history, summed once per step down in α.
# `sampled` takes flicker as a continuous 1/f process read every τ0 instead
# (structure function D(k) = 2·Δ²[t² ln t](k)).
function factor(alpha, n; sampled=false)
    F = Matrix{Float64}(I, n, n)
    if isodd(alpha) && sampled
        g(t) = t == 0 ? 0.0 : t^2 * log(abs(t))
        D(k) = 2 * (g(k + 1) - 2g(k) + g(k - 1))
        F = Matrix(cholesky(Symmetric([D(n) - D(abs(i - j)) / 2 for i in 1:n, j in 1:n])).L)
    elseif isodd(alpha)
        h = ones(4n); for k in 1:4n-1; h[k+1] = h[k] * (k - 0.5) / k; end
        H = [i >= j ? h[i-j+1] : 0.0 for i in 3n+1:4n, j in 1:4n]
        F = Matrix(cholesky(Symmetric(H * H')).L)
    end
    for _ in 1:(isodd(alpha) ? 1 - alpha : 2 - alpha)÷2; cumsum!(F, F, dims=1); end
    return F
end

# One window gives the bias. Flicker rows: the mean of the two flicker models.
function bias(a, m)
    d = vcat(ones(m), fill(-3.0, m), fill(3.0, m), fill(-1.0, m)) ./ sqrt(6.0 * m^4)
    B = [(F = factor(a, 4m; sampled); sum(abs2, window_op(m) * F) / sum(abs2, F' * d)) for sampled in (isodd(a) ? (false, true) : (false,))]
    return sum(B) / length(B), (maximum(B) / minimum(B) - 1) / 2
end
println("Bias B = E[MHTOTVAR] / E[MHVAR]\n\nnoise              m = 16     64    256   and at m = 2 … 15")
for (a, name) in NOISES
    println(rpad(name, 16), join((@sprintf("%7.4f", bias(a, m)[1]) for m in (16, 64, 256))),
            "   (", join((@sprintf("%.4f", bias(a, m)[1]) for m in 2:15), ", "), ")",
            isodd(a) ? @sprintf("   flicker models ±%.1f %% at m = 2", 100 * bias(a, 2)[2]) : "")
end

println("\nEDF = b·T/τ − c + d·τ/T + e·(τ/T)³, least squares on relative residuals, 16 ≤ m ≤ N/4, N = 512, 1024, 2048\n")
println("noise                 b       c       d       e   largest residual")
for (a, name) in NOISES
    s = Float64[]; y = Float64[]
    for N in (512, 1024, 2048)
        F = factor(a, N)
        for m in filter(<=(N ÷ 4), unique(round.(Int, 2 .^ (4:0.25:9))))
            L = 4m; n = N - L + 1
            q = (W = window_op(m); W' * W); Q = zeros(N, N)
            for o in 0:n-1; @views Q[o+1:o+L, o+1:o+L] .+= q; end
            G = F' * Q * F
            push!(s, m / (N - 1)); push!(y, tr(G)^2 / sum(abs2, G))
        end
    end
    X = hcat(1 ./ s, -ones(length(s)), s, s .^ 3)
    c = (X ./ y) \ ones(length(y))
    println(rpad(name, 16), @sprintf("%8.4f%8.3f%8.2f%8.1f%10.1f %%", c..., 100 * maximum(abs.(X * c ./ y .- 1))))
    a == 2 || continue
    # Below m = 16 the white PM model reads high: exact over model, mean of two record lengths.
    r = map((256, 1024)) do N
        map(2:15) do m
            n = N - 4m + 1; q = (W = window_op(m); W' * W); Q = zeros(N, N)
            for o in 0:n-1; @views Q[o+1:o+4m, o+1:o+4m] .+= q; end
            tr(Q)^2 / sum(abs2, Q) / ([(N - 1) / m, -1, m / (N - 1), (m / (N - 1))^3]' * c)
        end
    end
    println("  white PM, exact EDF over the model at m = 2 … 15: (", join((@sprintf("%.3f", v) for v in sum(r) ./ 2), ", "), ")")
end
