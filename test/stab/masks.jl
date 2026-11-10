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
    @test_throws ArgumentError telecom_mask(:G999, :mtie, 1.0)
    @test_throws ArgumentError telecom_mask(:G811, :adev, 1.0)

    # White PM of 1 ns rms: TDEV ≈ 1 ns/√(3m) passes the EEC mask; ×10 fails it.
    x = 1e-9 .* randn(Xoshiro(7), 30_000)
    c = mask_check(tdev(PhaseData(x, 1 / 30); ci=false), :G8262_opt1)
    @test c.pass && c.margin > 3 && length(c.limit) == length(c.tau)
    @test !mask_check(tdev(PhaseData(10 .* x, 1 / 30); ci=false), :G8262_opt1).pass
    @test mask_check(mtie(PhaseData(10 .* x, 1 / 30)), :G811).margin < 1   # ~80 ns p-p > 25 ns
    @test_throws ArgumentError mask_check(adev(PhaseData(x, 1.0)), :G811)
end
