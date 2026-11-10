# The PM bias and EDF models of TOTDEV and HTOTDEV against the exact moments
# of the kernels. A squared deviation is a quadratic form x'Qx; Q comes from
# the kernel by polarization on unit vectors, and for Gaussian noise with
# covariance Σ the mean is tr(QΣ) and the EDF tr(QΣ)² / tr((QΣ)²). The models
# are derived independently of the kernels in tools/total_pm_models.jl.

using SigmaTau: _totdev_core, _adev_core, _htotdev_core, _hdev_core, bias_correction, calculate_edf

@testset "TOTDEV and HTOTDEV on PM: bias and EDF are the kernels' exact moments" begin
    N = 128
    e(i) = Float64.((1:N) .== i)
    form(core, m) = (s = [core(e(i), [m], 1.0)[1]^2 for i in 1:N];
                     [(core(e(i) + e(j), [m], 1.0)[1]^2 - s[i] - s[j]) / 2 for i in 1:N, j in 1:N])
    # White PM, and flicker PM as a sampled 1/f process: generalized covariance
    # −D(k)/2 with structure function D(k) = 2·Δ²[t² ln t](k).
    g(t) = t == 0 ? 0.0 : t^2 * log(abs(t))
    Σ = Dict(:WHPM => Float64.((1:N) .== (1:N)'),
             :FLPM => [-(g(i - j + 1) - 2g(i - j) + g(i - j - 1)) for i in 1:N, j in 1:N])
    for (est, core, parent, var, ms) in ((:totdev, _totdev_core, _adev_core, :totvar, (16, 60)),
                                         (:htotdev, _htotdev_core, _hdev_core, :htot, (16, 40))),
        m in ms, noise in (:WHPM, :FLPM)
        A = form(core, m) * Σ[noise]
        mean, edf = sum(A[i, i] for i in 1:N), sum(A[i, i] for i in 1:N)^2 / sum(A .* A')
        B = mean / sum(form(parent, m) .* Σ[noise])
        @test bias_correction([noise], var, [m], N)[1] ≈ B rtol = 0.012
        @test calculate_edf(est, [1.0], [noise], [m], [Float64(m)], N, N - 1.0)[1] ≈ edf rtol = 0.08
    end
end
