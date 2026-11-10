# tools/mhtot_models.jl — exact bias and EDF of MHTOTDEV: the source of
# `_MHTOT_B` and `_coeff_mhtot` in src/edf.jl. No random numbers.
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
function factor(alpha, n)
    F = Matrix{Float64}(I, n, n)
    if isodd(alpha)
        h = ones(4n); for k in 1:4n-1; h[k+1] = h[k] * (k - 0.5) / k; end
        H = [i >= j ? h[i-j+1] : 0.0 for i in 3n+1:4n, j in 1:4n]
        F = Matrix(cholesky(Symmetric(H * H')).L)
    end
    for _ in 1:(isodd(alpha) ? 1 - alpha : 2 - alpha)÷2; cumsum!(F, F, dims=1); end
    return F
end

println("Bias B = E[MHTOTVAR] / E[MHVAR]\n\nnoise             m = 2      4      8     16     64    256")
for (a, name) in NOISES
    B = map((2, 4, 8, 16, 64, 256)) do m
        F = factor(a, 4m)
        d = vcat(ones(m), fill(-3.0, m), fill(3.0, m), fill(-1.0, m)) ./ sqrt(6.0 * m^4)
        sum(abs2, window_op(m) * F) / sum(abs2, F' * d)
    end
    println(rpad(name, 16), join((@sprintf("%7.4f", b) for b in B)))
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
end
