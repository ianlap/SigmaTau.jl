# test/stab/taus.jl — TauMode grid selector and tau_values.
# Loops collect the failing cases, so a failure names them.

using Test
using Random
using SigmaTau

@testset "TauMode + tau_values" begin
    kernels = (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
               :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
               :mtie, :pdev)
    modes = (AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade)

    @test isempty([(N, k) for N in (256, 1024, 4096), k in kernels
                   if tau_values(Octave, N, k) != SigmaTau._default_m_values(N, k)])
    @test isempty([(N, k) for N in (256, 1024), k in kernels
                   if tau_values(AllTaus, N, k) != 1:SigmaTau._kernel_m_max(N, k)])

    # Every grid is non-empty, sorted, unique, starts at 1 and stays in [1, m_max].
    valid(ms, mmax) = !isempty(ms) && issorted(ms) && allunique(ms) &&
                      first(ms) == 1 && all(m -> 1 <= m <= mmax, ms)
    @test isempty([(mode, N, k) for mode in modes, N in (128, 777, 4096), k in kernels
                   if !valid(tau_values(mode, N, k), SigmaTau._kernel_m_max(N, k))])

    @testset "grids at N = 1024 (:adev, m_max = 511)" begin
        @test SigmaTau._kernel_m_max(1024, :adev) == 511
        @test tau_values(Octave, 1024, :adev) == [1, 2, 4, 8, 16, 32, 64, 128, 256]
        @test tau_values(Decade, 1024, :adev) == [1, 10, 100]
        @test tau_values(HalfDecade, 1024, :adev) == [1, 3, 10, 32, 100, 316]
        @test first(tau_values(HalfOctave, 1024, :adev), 6) == [1, 2, 3, 4, 6, 8]
    end

    @testset "a TauMode argument equals the explicit tau_values form" begin
        pd = PhaseData(cumsum(randn(MersenneTwister(7), 1024)) .* 1e-9, 1.0)
        fns = (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
               htotdev, mhtotdev, mtie, pdev)
        @test isempty([k for (f, k) in zip(fns, kernels)
                       if f(pd, Decade; ci=false).dev !=
                          f(pd, tau_values(Decade, 1024, k); ci=false).dev])
        fd = FrequencyData(randn(MersenneTwister(8), 512) .* 1e-9, 1.0)
        @test adev(fd, Octave; ci=false).dev ==
              adev(fd, tau_values(Octave, 512, :adev); ci=false).dev
    end

    @test_throws ArgumentError tau_values(Octave, 1, :adev)     # N too short
    @test_throws ArgumentError tau_values(Octave, 1024, :nope)  # unknown kernel
end
