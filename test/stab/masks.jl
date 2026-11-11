@testset "telecom masks" begin
    # Breakpoints from the recommendations' tables, in ns.
    ns(m, s, τ) = telecom_mask(m, s, τ) * 1e9
    @test ns(:G811, :mtie, 1000.0) ≈ 300 && ns(:G811, :mtie, 2000.0) ≈ 310
    @test ns(:G811, :tdev, 1000.0) ≈ 30 && isnan(ns(:G811, :tdev, 2e4))
    @test ns(:G813_opt1, :mtie, 100.0) ≈ 40 * 100^0.1 && ns(:G813_opt1, :mtie, 1000.0) ≈ 25.25 * 1000^0.2
    @test ns(:G8262_opt1, :tdev, 25.0) ≈ 3.2 && ns(:G8262_opt1, :tdev, 100.0) ≈ 6.4
    @test ns(:G8262_opt2, :tdev, 1.0) ≈ 3.2 && ns(:G8262_opt2, :mtie, 10.0) ≈ 20 * 10^0.48
    @test ns(:G823_SSU, :mtie, 2000.0) ≈ 2000 && ns(:G823_SSU, :tdev, 100.0) ≈ 70
    @test ns(:G823_SEC, :tdev, 17.14) ≈ 12 && ns(:G823_PDH, :mtie, 20.0) ≈ 2000
    @test isnan(ns(:G8262_opt1, :mtie, 0.05))
    @test ns(:G811_1, :mtie, 1000.0) ≈ 15.0375 && ns(:G811_1, :mtie, 2000.0) ≈ 16.0375 && ns(:G811_1, :tdev, 1e4) ≈ 1
    @test ns(:G812_I, :mtie, 400.0) ≈ 160 && ns(:G812_I, :tdev, 50.0) ≈ 6
    @test ns(:G812_II_III, :mtie, 10.0) ≈ 40 * 10^0.4 && ns(:G812_II_III, :tdev, 1e5) ≈ 10
    @test ns(:G8262_1, :mtie, 0.1) ≈ 10 * 0.1^0.155 && ns(:G8262_1, :tdev, 50.0) ≈ 0.128 * √50     # 0.1 ≤ τ
    @test ns(:G8272_PRTC_A, :mtie, 273.0) ≈ 100.075 && ns(:G8272_PRTC_A, :mtie, 1e5) ≈ 100
    @test ns(:G8272_PRTC_B, :mtie, 54.5) ≈ 39.9875 && ns(:G8272_PRTC_B, :mtie, 1e5) ≈ 40
    @test ns(:G8272_PRTC_A, :tdev, 500.0) ≈ 15 && ns(:G8272_PRTC_B, :tdev, 500.0) ≈ 5
    @test ns(:G8272_1, :mtie, 100.0) ≈ 15.004 && ns(:G8272_1, :mtie, 4e5) ≈ 30 && ns(:G8272_1, :tdev, 6e4) ≈ 3.33333e-5 * 6e4
    # τ < 10 000 s in G.811 and for PRTC-A, τ ≤ 10 000 s at the G.823 PRC interface.
    @test isnan(ns(:G811, :tdev, 1e4)) && isnan(ns(:G8272_PRTC_A, :tdev, 1e4)) && ns(:G823_PRC, :tdev, 1e4) ≈ 30
    @test isnan(ns(:G8272_PRTC_B, :tdev, 1e5)) && isnan(ns(:G8272_1, :tdev, 1e6)) && ns(:G8272_1, :tdev, 999_999.0) ≈ 10
    @test_throws ArgumentError telecom_mask(:G999, :mtie, 1.0)
    @test_throws ArgumentError telecom_mask(:G811, :adev, 1.0)

    # White PM of 1 ns rms: TDEV ≈ 1 ns/√(3m) passes the EEC mask; ×10 fails it.
    x = 1e-9 .* randn(Xoshiro(7), 30_000)
    c = mask_check(tdev(PhaseData(x, 1 / 30); ci=false), :G8262_opt1)
    @test c.pass && c.margin > 3 && length(c.limit) == length(c.tau)
    @test !mask_check(tdev(PhaseData(10 .* x, 1 / 30); ci=false), :G8262_opt1).pass
    @test mask_check(mtie(PhaseData(10 .* x, 1 / 30)), :G811).margin < 1   # ~80 ns p-p > 25 ns
    @test mask_check(mtie(PhaseData(x, 1 / 30)), :G8272_PRTC_B).pass && !mask_check(mtie(PhaseData(10 .* x, 1 / 30)), :G8272_PRTC_A).pass
    @test_throws ArgumentError mask_check(adev(PhaseData(x, 1.0)), :G811)
end
