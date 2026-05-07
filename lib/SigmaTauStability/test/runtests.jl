using Test
using SigmaTauBase
using SigmaTauStability

@testset "SigmaTauStability Math Core Parity" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "ADEV Core" begin
        devs = _adev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
        # 2nd diff of x^2 is constant 2. 
        # For m=1, d2 = x[i+2] - 2x[i+1] + x[i] = 2. sum_sq = sum(4) = 4 * 98. 
        # dev = sqrt(4 * 98 / (2 * 98 * 1^2 * 1^2)) = sqrt(2) ≈ 1.414
        @test isapprox(devs[1], sqrt(2.0))
    end

    @testset "MDEV Core" begin
        devs = _mdev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
    end
    
    @testset "TDEV Core" begin
        devs = _tdev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
    end
    
    @testset "HDEV & MHDEV Core" begin
        # 3rd difference of x^2 is 0. So dev should be very close to 0.
        devs_h = _hdev_core(x, m_values, tau0)
        devs_mh = _mhdev_core(x, m_values, tau0)
        @test length(devs_h) == 3
        @test all(isfinite, devs_h)
        @test all(d -> d < 1e-10, devs_h)
        @test length(devs_mh) == 3
        @test all(isfinite, devs_mh)
        @test all(d -> d < 1e-10, devs_mh)
    end

    @testset "Total Deviations Core" begin
        # Just ensure they run and produce finite values.
        devs_tot = _totdev_core(x, m_values, tau0)
        devs_mtot = _mtotdev_core(x, m_values, tau0)
        devs_htot = _htotdev_core(x, m_values, tau0)
        devs_mhtot = _mhtotdev_core(x, m_values, tau0)
        
        @test length(devs_tot) == 3
        @test all(isfinite, devs_tot)
        @test length(devs_mtot) == 3
        @test all(isfinite, devs_mtot)
        @test length(devs_htot) == 3
        @test all(isfinite, devs_htot)
        @test length(devs_mhtot) == 3
        @test all(isfinite, devs_mhtot)
    end
    
    @testset "Typed API Wrappers & Stats (Lag1 and B1-ratio)" begin
        using Random
        Random.seed!(42)
        # N=100. For m=1, N_eff = 100 >= 50 -> Lag1 ACF used.
        # For m=4, N_eff = 25 < 50 -> B1-Ratio used.
        pd = PhaseData(x .+ randn(100), tau0)
        
        res_adev = adev(pd, m_values)
        @test length(res_adev.dev) == 3
        @test res_adev.noise_type[1] != :unknown
        @test all(isfinite, res_adev.ci_lower)
        @test all(isfinite, res_adev.ci_upper)
        
        res_hdev = hdev(pd, m_values)
        @test length(res_hdev.dev) == 3
        @test all(isfinite, res_hdev.ci_upper)
        
        res_tot = totdev(pd, m_values)
        @test length(res_tot.dev) == 3
        @test all(isfinite, res_tot.ci_upper)
        
        res_ldev = ldev(pd, m_values)
        @test length(res_ldev.dev) == 3
    end
end
