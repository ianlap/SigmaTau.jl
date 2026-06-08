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
        @test r2.ci_lower  ≈ r.ci_lower
        @test r2.ci_upper  ≈ r.ci_upper
        @test r2.edf       ≈ r.edf

        # Round-trip without CI — empty vectors must survive the cycle
        r_nci = adev(pd, ms; ci=false)
        save_result(path, r_nci)
        r3 = load_result(path)
        @test r3.deviation_type === :adev
        @test r3.tau ≈ r_nci.tau
        @test r3.dev ≈ r_nci.dev
        @test isempty(r3.noise_type)
        @test isempty(r3.ci_lower)
        @test isempty(r3.ci_upper)
        @test isempty(r3.edf)

        @test save_result(path, r) == path
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
            @test s2[k].ci_lower   ≈ s[k].ci_lower
            @test s2[k].ci_upper   ≈ s[k].ci_upper
            @test s2[k].edf        ≈ s[k].edf
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
        @test isempty(s2[:adev].ci_lower)
        @test isempty(s2[:mtie].ci_lower)
        @test s2[:adev].dev ≈ s[:adev].dev
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
