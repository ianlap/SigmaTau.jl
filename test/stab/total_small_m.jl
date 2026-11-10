# The bias factors of the modified and Hadamard totals below m = 16, and the
# white PM EDF of MHTOTDEV there, against the exact moments of the kernels.
# For x = F·w with white w, E[σ²] is the kernel's squared output summed over
# the columns of F; on white PM, F = I, the EDF is tr(Q)² / ‖Q‖² with Q the
# kernel's form by polarization. Derived independently of the kernels in
# tools/total_pm_models.jl and tools/mhtot_models.jl.

using SigmaTau: _mtotdev_core, _mdev_core, _htotdev_core, _hdev_core, _mhtotdev_core, _mhdev_core,
                bias_correction, calculate_edf

@testset "Total-family bias and EDF below m = 16 are the kernels' exact moments" begin
    N = 128
    white = [Float64.((1:N) .== j) for j in 1:N]    # white PM
    walk = cumsum.(white)                           # white FM
    sq(core, x, m) = core(x, [m], 1.0)[1]^2
    for (var, core, parent, noise, F) in ((:mtot, _mtotdev_core, _mdev_core, :WHPM, white),
                                          (:htot, _htotdev_core, _hdev_core, :WHFM, walk),
                                          (:mhtot, _mhtotdev_core, _mhdev_core, :WHPM, white)), m in (2, 3)
        @test bias_correction([noise], var, [m], N)[1] ≈ sum(sq.(core, F, m)) / sum(sq.(parent, F, m)) rtol = 2e-4
    end

    s = sq.(_mhtotdev_core, white, 2)
    Q = [(sq(_mhtotdev_core, white[i] + white[j], 2) - s[i] - s[j]) / 2 for i in 1:N, j in 1:N]
    @test calculate_edf(:mhtotdev, [1.0], [:WHPM], [2], [2.0], N, N - 1.0)[1] ≈ sum(s)^2 / sum(abs2, Q) rtol = 0.03
end
