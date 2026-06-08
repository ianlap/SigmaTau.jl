# test/stab/taus.jl — TauMode grid selector + tau_values.
# Standalone-runnable (`julia --project=. test/stab/taus.jl`) and included
# from test/stab/runtests.jl.

using Test
using Random
using SigmaTau

@testset "TauMode + tau_values" begin
    kernels = (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
               :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
               :mtie, :pdev)

    @testset "Octave matches _default_m_values (byte-identical default)" begin
        for N in (256, 1024, 4096), k in kernels
            @test tau_values(Octave, N, k) == SigmaTau._default_m_values(N, k)
        end
    end

    @testset "AllTaus is 1:m_max" begin
        for N in (256, 1024), k in kernels
            mmax = SigmaTau._kernel_m_max(N, k)
            @test tau_values(AllTaus, N, k) == collect(1:mmax)
        end
    end

    @testset "Specific grids at N=1024 (:adev, m_max=511)" begin
        @test SigmaTau._kernel_m_max(1024, :adev)        == 511
        @test tau_values(Octave, 1024, :adev)            == [1, 2, 4, 8, 16, 32, 64, 128, 256]
        @test tau_values(Decade, 1024, :adev)            == [1, 10, 100]
        @test tau_values(HalfDecade, 1024, :adev)        == [1, 3, 10, 32, 100, 316]
        @test first(tau_values(HalfOctave, 1024, :adev), 6) == [1, 2, 3, 4, 6, 8]
    end

    @testset "All modes: non-empty, sorted, unique, within [1, m_max]" begin
        for mode in (AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade)
            for N in (128, 777, 4096), k in kernels
                ms   = tau_values(mode, N, k)
                mmax = SigmaTau._kernel_m_max(N, k)
                @test !isempty(ms)
                @test issorted(ms)
                @test allunique(ms)
                @test first(ms) == 1
                @test all(m -> 1 <= m <= mmax, ms)
            end
        end
    end

    @testset "Deviation TauMode method == explicit tau_values form" begin
        Random.seed!(7)
        pd = PhaseData(cumsum(randn(1024)) .* 1e-9, 1.0)
        for (f, k) in ((adev, :adev), (mdev, :mdev), (tdev, :tdev), (hdev, :hdev),
                       (mhdev, :mhdev), (htdev, :htdev), (totdev, :totdev),
                       (mtotdev, :mtotdev), (ttotdev, :ttotdev), (htotdev, :htotdev),
                       (mhtotdev, :mhtotdev), (mtie, :mtie), (pdev, :pdev))
            a = f(pd, Decade; ci=false)
            b = f(pd, tau_values(Decade, length(pd.x), k); ci=false)
            @test a.tau == b.tau
            @test a.dev == b.dev
        end
        # FrequencyData TauMode path delegates through _freq_to_phase.
        fd = FrequencyData(randn(512) .* 1e-9, 1.0)
        @test adev(fd, Octave; ci=false).dev ==
              adev(fd, tau_values(Octave, length(fd.y), :adev); ci=false).dev
    end

    @testset "error paths" begin
        @test_throws ArgumentError tau_values(Octave, 1, :adev)   # N too short
        @test_throws ArgumentError tau_values(Octave, 1024, :nope) # unknown kernel
    end
end
