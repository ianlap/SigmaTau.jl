# test/io/results.jl — save_result/load_result (v1) and save_suite/load_suite (v2).

@testset "Result & suite I/O" begin
    tmpdir = mktempdir()

    @testset "StabilityResult v1 round-trip" begin
        path = joinpath(tmpdir, "result.tsv")
        Random.seed!(20260510)
        pd = PhaseData(randn(128) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8]

        # Round-trip with CI
        r  = adev(pd, ms; ci=true)
        save_result(path, r)
        r2 = load_result(path)
        @test r2.deviation_type === r.deviation_type
        @test r2.tau       ≈ r.tau
        @test r2.dev       ≈ r.dev
        @test r2.noise_type == r.noise_type
        @test ci_lower(r2) ≈ ci_lower(r)
        @test ci_upper(r2) ≈ ci_upper(r)
        @test r2.edf       ≈ r.edf
        @test r2.neff      == r.neff
        @test r2.neff      == 128 .- 2 .* ms    # adev: N − 2m windows

        # Round-trip without CI — empty vectors must survive the cycle,
        # but neff is always present.
        r_nci = adev(pd, ms; ci=false)
        save_result(path, r_nci)
        r3 = load_result(path)
        @test r3.deviation_type === :adev
        @test r3.tau ≈ r_nci.tau
        @test r3.dev ≈ r_nci.dev
        @test isempty(r3.noise_type)
        @test isempty(r3.ci)
        @test isempty(r3.edf)
        @test r3.neff == r_nci.neff

        @test save_result(path, r) == path
    end

    @testset "legacy 6-column file (no neff) still loads" begin
        # Files written before the neff column existed carry six data columns;
        # they load with neff filled with zeros and everything else intact.
        path = joinpath(tmpdir, "legacy.tsv")
        open(path, "w") do io
            println(io, "# SigmaTau StabilityResult v1")
            println(io, "# deviation_type=adev")
            println(io, "# calc_ci=true")
            println(io, "tau\tdev\tnoise_type\tci_lower\tci_upper\tedf")
            println(io, "1.0\t1.5e-9\tWHFM\t1.4e-9\t1.7e-9\t98.5")
            println(io, "2.0\t1.1e-9\tWHFM\t1.0e-9\t1.3e-9\t60.2")
        end
        r = load_result(path)
        @test r.deviation_type === :adev
        @test r.tau == [1.0, 2.0]
        @test r.dev == [1.5e-9, 1.1e-9]
        @test r.ci[1] == (lo=1.4e-9, hi=1.7e-9)
        @test ci_upper(r) == [1.7e-9, 1.3e-9]
        @test r.edf == [98.5, 60.2]
        @test r.neff == [0, 0]
    end

    @testset "StabilitySuite v2 round-trip (with CI + metadata)" begin
        path = joinpath(tmpdir, "suite.tsv")
        Random.seed!(20260511)
        pd = PhaseData(cumsum(randn(1024)) .* 1e-9, 1.0)
        s  = stability(pd; devs=(:adev, :mdev, :hdev), taus=Octave, ci=true, confidence=0.9)

        @test save_suite(path, s; source_file="clock_a.dat") == path
        s2 = load_suite(path)
        @test s2 isa StabilitySuite
        @test keys(s2) == keys(s)
        @test s2.data_kind === s.data_kind
        @test s2.tau0 == s.tau0
        @test s2.n == s.n
        @test s2.confidence == s.confidence
        @test s2.tau_mode === s.tau_mode
        for k in keys(s)
            @test s2[k].tau        ≈ s[k].tau
            @test s2[k].dev        ≈ s[k].dev
            @test s2[k].noise_type == s[k].noise_type
            @test ci_lower(s2[k])  ≈ ci_lower(s[k])
            @test ci_upper(s2[k])  ≈ ci_upper(s[k])
            @test s2[k].edf        ≈ s[k].edf
            @test s2[k].neff       == s[k].neff
        end

        # Metadata header is present and parseable.
        hdr = readlines(path)
        @test occursin("StabilitySuite v2", hdr[1])
        @test any(l -> startswith(l, "# source_file=clock_a.dat"), hdr)
        @test any(l -> startswith(l, "# package_version="), hdr)
        @test any(l -> startswith(l, "# timestamp="), hdr)
    end

    @testset "suite round-trip without CI" begin
        path = joinpath(tmpdir, "suite_nci.tsv")
        Random.seed!(20260512)
        pd = PhaseData(cumsum(randn(256)) .* 1e-9, 1.0)
        s  = stability(pd; devs=(:adev, :mtie), taus=Octave, ci=false)
        save_suite(path, s)
        s2 = load_suite(path)
        @test isnan(s2.confidence)
        @test isempty(s2[:adev].ci)
        @test isempty(s2[:mtie].ci)
        @test s2[:adev].dev ≈ s[:adev].dev
        @test s2[:adev].neff == s[:adev].neff
        @test s2[:mtie].neff == s[:mtie].neff
    end

    @testset "cross-format guards" begin
        rpath = joinpath(tmpdir, "single.tsv")
        spath = joinpath(tmpdir, "multi.tsv")
        Random.seed!(20260513)
        pd = PhaseData(cumsum(randn(128)) .* 1e-9, 1.0)
        save_result(rpath, adev(pd, [1, 2, 4]; ci=false))
        save_suite(spath, stability(pd; devs=(:adev,), ci=false))
        @test_throws ErrorException load_suite(rpath)    # v1 file via load_suite
        @test_throws ErrorException load_result(spath)   # v2 file via load_result
    end
end
