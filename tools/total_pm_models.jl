# tools/total_pm_models.jl — exact bias and EDF of TOTDEV and HTOTDEV on white
# and flicker PM: the source of the PM rows of `_total_bias`, `_totvar_pm_edf`
# and `_htot_pm_edf` in src/edf.jl. No random numbers.
#
#   julia -t auto tools/total_pm_models.jl        # under a minute
#
# Both statistics are quadratic forms x'Qx in the phase record, so for
# Gaussian noise with covariance Σ = FF' the mean is tr(F'QF) and
# EDF = tr(F'QF)² / ‖F'QF‖². The forms are built from the definitions,
# independently of the package kernels; test/stab/total_pm_models.jl holds
# the kernels to the same numbers. Linear coefficients are fitted here; the
# other constants of the EDF forms (CONST) were set by eye and their
# residuals are printed.
using LinearAlgebra, Printf

const CONST = (kappa=0.93, c_flpm=4.75,
               htot=Dict(2 => (c=2.0, s0=0.125, a=1.0, b=0.03), 1 => (c=2.5, s0=0.16, a=0.6, b=0.04)))

# F with Σ = FF' for n phase samples. White PM is the identity. Flicker PM is
# either a continuous 1/f process sampled at τ0 (structure function
# D(k) = 2·Δ²[t² ln t](k), the model of the Greenhall–Riley EDF algorithm) or
# Kasdin's half-order filter over a 4n history (the model of `noise_gen`).
function factor(noise, n)
    noise == :white && return Matrix{Float64}(I, n, n)
    if noise == :sampled
        g(t) = t == 0 ? 0.0 : t^2 * log(abs(t))
        D(k) = 2 * (g(k + 1) - 2g(k) + g(k - 1))
        return Matrix(cholesky(Symmetric([D(n) - D(abs(i - j)) / 2 for i in 1:n, j in 1:n])).L)
    end
    h = ones(4n); for k in 1:4n-1; h[k+1] = h[k] * (k - 0.5) / k; end
    H = [i >= j ? h[i-j+1] : 0.0 for i in 3n+1:4n, j in 1:4n]
    return Matrix(cholesky(Symmetric(H * H')).L)
end

# TOTVAR = ‖Mx‖²: second differences at lag m of the record extended by odd
# reflection about each end sample.
function tot_op(N, m)
    M = zeros(N - 2, N)
    put!(r, k, c) = k < 1 ? (r[1] += 2c; r[2-k] -= c) : k > N ? (r[N] += 2c; r[2N-k] -= c) : (r[k] += c)
    for n in 2:N-1
        r = @view M[n-1, :]
        put!(r, n - m, 1.0); put!(r, n, -2.0); put!(r, n + m, 1.0)
    end
    return M ./ sqrt(2.0 * (N - 2) * m^2)
end

# One TOTHVAR window, 3m frequency samples from 3m + 1 phase samples: remove
# the half-average slope, mirror at both ends, and take the second difference
# of m-sample sums at 6m positions.
function htot_window(m)
    L = 3m
    D = zeros(L, L + 1); for j in 1:L; D[j, j] = -1; D[j, j+1] = 1; end
    h = L ÷ 2; lo = cld(L, 2) + 1
    slope = zeros(L); slope[1:h] .= -1 / h; slope[lo:L] .= 1 / (L - lo + 1)
    slope ./= isodd(L) ? (L + 1) / 2 : L / 2
    P = I - ((0:L-1) .- Float64(h)) * slope'
    E = zeros(3L, L); for k in 1:L; E[k, L+1-k] = 1; E[L+k, k] = 1; E[2L+k, L+1-k] = 1; end
    B = zeros(6m, 3L)
    for j in 1:6m; B[j, j:j+m-1] .= 1; B[j, j+m:j+2m-1] .= -2; B[j, j+2m:j+3m-1] .= 1; end
    return B * E * P * D ./ sqrt(36.0 * m^3)
end

function htot_G(F, m)
    N = size(F, 1); W = htot_window(m); q = W' * W; n = N - 3m; Q = zeros(N, N)
    for o in 0:n-1; @views Q[o+1:o+3m+1, o+1:o+3m+1] .+= q; end
    return F' * Q * F ./ n
end

# Expectation of the parent on the same noise: AVAR (d = 2) or HVAR (d = 3).
function parent(noise, m, d)
    c = d == 2 ? [1.0, -2, 1] : [-1.0, 3, -3, 1]
    v = zeros(d * m + 1); v[1:m:end] .= c
    return sum(abs2, factor(noise, d * m + 1)' * v) / (sum(abs2, c) / (d == 2 ? 3 : 10 / 3) * m^2)
end

edf(G) = tr(G)^2 / sum(abs2, G)
lsq(X, y) = (X ./ y) \ ones(length(y))                 # least squares on the relative residual
worst(model, exact) = 100 * maximum(abs.(model ./ exact .- 1))
grid(lo, hi) = filter(m -> lo <= m <= hi, unique(vcat(2:16, round.(Int, 2 .^ (4:0.25:12)))))

Ns = (512, 1024)
F = Dict((noise, N) => factor(noise, N) for noise in (:white, :sampled), N in Ns)

# ── TOTVAR ────────────────────────────────────────────────────────────────────
tot = Dict((noise, N, m) => (A = tot_op(N, m) * F[(noise, N)]; G = A * A'; (tr(G) / parent(noise, m, 2), edf(G)))
           for noise in (:white, :sampled), N in Ns for m in grid(2, N ÷ 2))
println("TOTVAR, white PM")
k = [key for key in keys(tot) if key[1] == :white]
@printf("  B = 1 + (8/3)·⌊m/2⌋/(N − 2): largest difference %.1e\n",
        maximum(abs(tot[(n, N, m)][1] - 1 - (8 / 3) * (m ÷ 2) / (N - 2)) for (n, N, m) in k))
@printf("  EDF = 72N²B²/(140N + 225m + 64m²): largest residual %.1f %%\n",
        worst([72.0 * N^2 * tot[(n, N, m)][1]^2 / (140N + 225m + 64.0 * m^2) for (n, N, m) in k], [tot[key][2] for key in k]))

# B²/EDF = g(m)/N + κ·(m/N)² at two record lengths gives g and κ for each m.
ms = grid(3, 256)
gk = [[1/512 (m/512)^2; 1/1024 (m/1024)^2] \ [tot[(:sampled, N, m)][1]^2 / tot[(:sampled, N, m)][2] for N in Ns] for m in ms]
L = log.(ms)
cg = lsq(hcat(ones(length(ms)), L, L .^ 2), ms ./ first.(gk))
g(m) = m / (cg[1] + cg[2] * log(m) + cg[3] * log(m)^2)
println("\nFlicker PM (sampled process)")
@printf("  AVAR interior, N/EDF → g(m) = m/(%.3f + %.3f ln m + %.3f ln² m): largest residual %.1f %%\n", cg..., worst(g.(ms), first.(gk)))
@printf("  κ for m ≥ 16: %.3f to %.3f\n", extrema(last.(gk)[ms .>= 16])...)
k = [key for key in keys(tot) if key[1] == :sampled]
@printf("  TOTVAR EDF = B²/(g/N + %.2f·(m/N)²): largest residual %.1f %%\n", CONST.kappa,
        worst([tot[(n, N, m)][1]^2 / (g(m) / N + CONST.kappa * (m / N)^2) for (n, N, m) in k], [tot[key][2] for key in k]))

# TOTVAR bias on flicker PM, the mean of the two models: B = 1 + k(m)·m/(N − 1).
ms = grid(8, 512)
Fk = factor(:kasdin, 1024)
kk(Fn, noise, m) = (sum(abs2, tot_op(1024, m) * Fn) / parent(noise, m, 2) - 1) * 1023 / m
ks = [kk(F[(:sampled, 1024)], :sampled, m) for m in ms]; kd = [kk(Fk, :kasdin, m) for m in ms]
ck = hcat(ones(length(ms)), log.(ms)) \ ((ks .+ kd) ./ 2)
@printf("  TOTVAR B = 1 + (%.3f + %.3f ln m)·τ/T: largest error in B %.1f %% (sampled), %.1f %% (Kasdin)\n", ck...,
        (worst(1 .+ (ck[1] .+ ck[2] .* log.(ms)) .* ms ./ 1023, 1 .+ kx .* ms ./ 1023) for kx in (ks, kd))...)

# ── TOTHVAR ───────────────────────────────────────────────────────────────────
# The bias depends on m alone: one window is enough.
ms = grid(2, 512)
hb(noise, m) = (W = htot_window(m); sum(abs2, W * factor(noise, 3m + 1)) / parent(noise, m, 3))
Bw = hb.(:white, ms); Bf = (hb.(:sampled, ms) .+ hb.(:kasdin, ms)) ./ 2
odd = isodd.(ms)
cw = lsq(hcat(odd, -odd ./ ms, .!odd, -.!odd ./ ms), Bw)
println("\nTOTHVAR bias")
@printf("  white PM: B = %.4f − %.4f/m (odd m), %.4f − %.4f/m (even m): largest residual %.2f %%\n", cw...,
        worst(hcat(odd, -odd ./ ms, .!odd, -.!odd ./ ms) * cw, Bw))
Lc = CONST.c_flpm .+ log.(ms)
X = hcat(1 ./ Lc, log.(ms) ./ Lc, -odd ./ ms, -.!odd ./ ms)
cf = lsq(X, Bf)
@printf("  flicker PM: B = (%.3f + %.3f ln m)/(%.2f + ln m) − %.3f/m (odd m) or %.3f/m (even m)\n", cf[1], cf[2], CONST.c_flpm, cf[3], cf[4])
@printf("    largest residual, m ≥ 16: %.1f %% (sampled), %.1f %% (Kasdin); m = 2: %+.1f %%, %+.1f %%\n",
        (worst((X*cf)[ms .>= 16], hb.(n, ms[ms .>= 16])) for n in (:sampled, :kasdin))..., (100 * ((X*cf)[1] / hb(n, 2) - 1) for n in (:sampled, :kasdin))...)

# EDF = B²·(N/g)·(1 − a·exp(−(N − 3m)/(bN))) / (1 + c·min(s, s0)), s = m/(N − 1),
# with g = 2.45 on white PM and the AVAR interior g(m) on flicker PM.
println("\nTOTHVAR EDF, residual of the form by τ/T (N = 1024), and its largest for m ≥ 16 at both lengths")
for (alpha, noise, gg) in ((2, :white, m -> 2.45), (1, :sampled, g))
    p = CONST.htot[alpha]
    form(N, m, B) = B^2 * (N / gg(m)) * (1 - p.a * exp(-(N - 3m) / (p.b * N))) / (1 + p.c * min(m / (N - 1), p.s0))
    res = Dict((N, m) => (G = htot_G(F[(noise, N)], m); form(N, m, tr(G) / parent(noise, m, 3)) / edf(G) - 1)
               for N in Ns for m in unique(vcat(grid(2, (N - 1) ÷ 3), round.(Int, range(0.25N, (N - 1) / 3, length=12)))))
    near(s) = argmin(m -> abs(m / 1023 - s), [m for (N, m) in keys(res) if N == 1024])
    println("  ", rpad(noise == :white ? "white PM" : "flicker PM", 12),
            join((@sprintf("%.3f: %+.0f %%", m / 1023, 100 * res[(1024, m)]) for m in near.((1 / 64, 1 / 16, 1 / 8, 1 / 4, 0.29, 0.31, 0.325, 1 / 3))), "   "),
            @sprintf("   largest %.0f %%", 100 * maximum(abs(v) for ((N, m), v) in res if m >= 16)))
end
