# test/io/results.jl — save/load of results (v1 format) and suites (v2).

# Every stored field survives the round trip.
roundtrips(a, b) =
    a.deviation_type === b.deviation_type && a.noise_type == b.noise_type &&
    a.neff == b.neff &&
    all(isapprox(p, q; nans=true) for (p, q) in ((a.tau, b.tau), (a.dev, b.dev),
        (ci_lower(a), ci_lower(b)), (ci_upper(a), ci_upper(b)), (a.edf, b.edf)))

@testset "Result & suite I/O" begin
    tmpdir = mktempdir()
    pd = PhaseData(cumsum(randn(MersenneTwister(20260510), 1024)) .* 1e-9, 1.0)

    @testset "StabilityResult round-trip, with and without CI" begin
        path = joinpath(tmpdir, "result.tsv")
        for ci in (true, false)
            r = adev(pd, [1, 2, 4, 8]; ci)
            @test save(path, r) == path
            @test roundtrips(load(path), r)
        end
        @test isempty(load(path).ci)
    end

    @testset "v0.4 file without the neff column loads with neff = 0" begin
        path = joinpath(tmpdir, "legacy.tsv")
        write(path, """
            # SigmaTau StabilityResult v1
            # deviation_type=adev
            # calc_ci=true
            tau\tdev\tnoise_type\tci_lower\tci_upper\tedf
            1.0\t1.5e-9\tWHFM\t1.4e-9\t1.7e-9\t98.5
            2.0\t1.1e-9\tWHFM\t1.0e-9\t1.3e-9\t60.2
            """)
        r = load(path)
        @test (r.deviation_type, r.tau, r.dev, r.ci[1], ci_upper(r), r.edf, r.neff) ==
              (:adev, [1.0, 2.0], [1.5e-9, 1.1e-9], (lo=1.4e-9, hi=1.7e-9),
               [1.7e-9, 1.3e-9], [98.5, 60.2], [0, 0])
    end

    @testset "StabilitySuite round-trip, with and without CI" begin
        path = joinpath(tmpdir, "suite.tsv")
        for ci in (true, false)
            s = stability(pd; devs=(:adev, :mdev, :mtie), m=Octave, ci, confidence=0.9)
            @test save(path, s; source_file="clock_a.dat") == path
            s2 = load(path)
            @test s2 isa StabilitySuite
            @test (keys(s2), s2.data_kind, s2.tau0, s2.n, s2.tau_mode) ==
                  (keys(s), s.data_kind, s.tau0, s.n, s.tau_mode)
            @test isequal(s2.confidence, s.confidence)
            @test all(roundtrips(s2[k], s[k]) for k in keys(s))
        end
        hdr = readlines(path)
        @test occursin("StabilitySuite v2", hdr[1])
        @test all(p -> any(startswith(p), hdr),
                  ("# source_file=clock_a.dat", "# package_version=", "# timestamp="))
    end

    @testset "load picks the file type" begin
        rpath, spath = joinpath(tmpdir, "single.tsv"), joinpath(tmpdir, "multi.tsv")
        save(rpath, adev(pd, [1, 2, 4]; ci=false))
        save(spath, stability(pd; devs=(:adev,), ci=false))
        @test load(rpath) isa StabilityResult
        @test load(spath) isa StabilitySuite
        @test_throws ErrorException SigmaTau._load_suite(rpath)
        @test_throws ErrorException SigmaTau._load_result(spath)
        save(rpath, pd)
        @test load(rpath).x == pd.x
        save(rpath, FrequencyData(diff(pd.x), 2.0))
        @test load(rpath) isa FrequencyData && load(rpath).tau0 == 2.0
    end
end
