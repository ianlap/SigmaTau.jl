# The shipped MHTOTDEV bias and EDF models against the exact moments of the
# kernel. The squared deviation is a quadratic form x'Qx, so for Gaussian
# noise x = F·w with white w, its mean is tr(G) and its EDF tr(G)² / ‖G‖²
# with G = F'QF, which the kernel itself yields by polarization. The models
# are derived independently of the kernel in tools/mhtot_models.jl.

using SigmaTau: _mhtotdev_core, _mhdev_core, bias_correction, calculate_edf

@testset "MHTOTDEV bias and EDF models are the kernel's exact moments" begin
    N = 129    # τ = 32 τ0 is T/4, the reach limit
    sq(core, x, m) = core(x, [m], 1.0)[1]^2
    for (noise, sums) in ((:WHFM, 1), (:RRFM, 3)), m in (16, 32)
        # Columns: the phase response to one white sample, summed `sums` times.
        F = collect.(eachcol(foldl((A, _) -> cumsum(A, dims=1), 1:sums; init=Float64.((1:N) .== (1:N)'))))
        s = sq.(_mhtotdev_core, F, m)
        G = [(sq(_mhtotdev_core, F[i] + F[j], m) - s[i] - s[j]) / 2 for i in 1:N, j in 1:N]

        B = sum(s) / sum(sq.(_mhdev_core, F, m))
        @test B ≈ bias_correction([noise], :mhtot, [m], N)[1] rtol = 0.004
        edf = calculate_edf(:mhtotdev, [1.0], [noise], [m], [Float64(m)], N, N - 1.0)[1]
        @test edf ≈ sum(s)^2 / sum(abs2, G) rtol = 0.03
    end
end
