# Tests for the sigmatau command line, run in process:
#
#   julia --project=cli cli/test/runtests.jl
#
# With SIGMATAU_EXE set to a compiled binary (cli/build/bin/sigmatau), the
# last testset also runs it and checks that it answers exactly like the
# in-process call.

using Test
using SigmaTau
using SigmaTauCLI

function run_cli(args...)
    o = IOBuffer()
    e = IOBuffer()
    code = SigmaTauCLI.sigmatau(collect(String, args); out=o, err=e)
    return (code=code, out=String(take!(o)), err=String(take!(e)))
end

# Deterministic phase record (seconds), long enough for octave grids.
const X = cumsum(sinpi.((1:256) ./ 7.3)) ./ 1e9

struct BrokenPipeIO <: IO end
Base.write(::BrokenPipeIO, ::UInt8) = throw(Base.IOError("write: broken pipe (EPIPE)", Base.UV_EPIPE))
Base.unsafe_write(::BrokenPipeIO, ::Ptr{UInt8}, ::UInt) = throw(Base.IOError("write: broken pipe (EPIPE)", Base.UV_EPIPE))

function write_table(path, x; tau0=1.0)
    open(path, "w") do io
        for (i, v) in enumerate(x)
            println(io, (i - 1) * tau0, ",", v)
        end
    end
    return path
end

csv_rows(s) = [split(l, ',') for l in split(strip(s), '\n')[2:end]]

@testset "sigmatau CLI" begin

@testset "version and help" begin
    for v in (["version"], ["--version"], ["-V"])
        r = run_cli(v...)
        @test r.code == 0
        @test occursin(string(pkgversion(SigmaTau)), r.out)
    end
    for h in (["help"], ["--help"], ["-h"], String[])
        r = run_cli(h...)
        @test r.code == 0
        for c in ("analyze", "check", "spectrum", "dynamic", "completion")
            @test occursin(c, r.out)
        end
    end
    for c in ("analyze", "check", "spectrum", "dynamic")
        r = run_cli(c, "--help")
        @test r.code == 0
        @test occursin("sigmatau " * c, r.out)
        @test occursin("--format", r.out)
        @test run_cli("help", c).out == r.out
    end
    @test occursin("--devs", run_cli("analyze", "--help").out)
    @test !occursin("--devs", run_cli("spectrum", "--help").out)
    @test occursin("--window", run_cli("dynamic", "--help").out)
end

@testset "unknown command or option" begin
    r = run_cli("frobnicate")
    @test r.code == 1
    @test occursin("unknown command", r.err)
    @test run_cli("--bogus").code == 1
    r = run_cli("analyze", "--synth", "64", "--bogus")
    @test r.code == 1
    @test occursin("unknown option --bogus", r.err)
    r = run_cli("check", "--synth", "64", "--devs", "adev")
    @test r.code == 1
    @test occursin("does not apply", r.err)
end

@testset "completion" begin
    for sh in ("bash", "zsh", "fish", "powershell")
        r = run_cli("completion", sh)
        @test r.code == 0
        s = r.out
        @test !occursin("@@", s)
        optpat = sh == "fish" ? "-l " : "--"
        for c in ("analyze", "check", "spectrum", "dynamic")
            @test occursin(c, s)
        end
        @test occursin(optpat * "devs", s)
        @test occursin(optpat * "out", s)
        @test occursin("mhtotdev", s) && occursin("dhdev", s)
        @test occursin("half-octave", s) && occursin("json", s)
    end
    @test occursin("-o default", run_cli("completion", "bash").out)
    @test occursin("_files", run_cli("completion", "zsh").out)
    @test occursin("-r -F", run_cli("completion", "fish").out)
    @test occursin("Register-ArgumentCompleter", run_cli("completion", "powershell").out)
    r = run_cli("completion", "tcsh")
    @test r.code == 1
    @test occursin("tcsh", r.err)
    @test run_cli("completion").code == 1
end

mktempdir() do dir
    pf = write_table(joinpath(dir, "phase.csv"), X)

    @testset "analyze -> CSV" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev,mdev")
        @test r.code == 0
        lines = split(strip(r.out), '\n')
        @test lines[1] == "source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type"
        rows = csv_rows(r.out)
        @test all(f -> length(f) == 10, rows)
        @test any(f -> f[3] == "adev", rows) && any(f -> f[3] == "mdev", rows)
        @test all(f -> f[2] == "phase", rows)
    end

    @testset "numbers match the library" begin
        pd = PhaseData(X, 1.0)
        ref = adev(pd, tau_values(Octave, length(X), :adev))
        rows = csv_rows(run_cli("analyze", pf, "--tau0", "1", "--devs", "adev").out)
        @test parse.(Float64, getindex.(rows, 4)) == ref.tau
        @test parse.(Float64, getindex.(rows, 5)) == ref.dev
        @test parse.(Float64, getindex.(rows, 6)) == [c.lo for c in ref.ci]
        @test getindex.(rows, 10) == String.(ref.noise_type)
        ref = totdev(pd, [1, 3, 9]; correct_bias=false)
        rows = csv_rows(run_cli("analyze", pf, "--devs", "totdev", "--m", "9,1,3", "--no-bias-correction").out)
        @test parse.(Float64, getindex.(rows, 5)) == ref.dev
    end

    @testset "frequency input" begin
        ff = write_table(joinpath(dir, "freq.csv"), diff(X))
        r = run_cli("analyze", ff, "--freq", "--tau0", "1", "--devs", "adev")
        @test r.code == 0
        @test all(f -> f[2] == "frequency", csv_rows(r.out))
    end

    @testset "analyze -> JSON" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev,tdev", "--format", "json")
        @test r.code == 0
        for k in ("\"schema_version\":1", "\"command\":\"analyze\"", "\"analysis\"", "\"results\"",
                  "\"points\"", "\"deviations\"", "\"noise_type\"", "\"tau_mode\":\"octave\"",
                  "\"axis\":\"sigma_x\"", "\"label\":\"TDEV\"", "\"record\"", "\"elapsed_ms\"")
            @test occursin(k, r.out)
        end
        @test !occursin("\"spectrum\"", r.out)
        @test !occursin("\"series\"", r.out)
        r = run_cli("analyze", pf, "--format=json", "--spectrum", "--series", "10")
        @test r.code == 0
        @test occursin("\"spectrum\":{\"nperseg\":64", r.out)
        @test occursin("\"series\":{\"t\":[0.0,26.0,", r.out)
    end

    @testset "no CI: nulls in JSON, blanks in CSV, never NaN" begin
        r = run_cli("analyze", pf, "--devs", "adev", "--format", "json", "--no-ci")
        @test r.code == 0
        @test occursin("\"ci\":false,\"confidence\":null", r.out)
        @test occursin("\"ci_lower\":null", r.out)
        @test !occursin("NaN", r.out)
        r = run_cli("analyze", pf, "--devs", "adev", "--no-ci")
        @test !occursin("NaN", r.out)
        for f in csv_rows(r.out)
            @test f[6] == "" && f[7] == "" && f[8] == "" && f[10] == ""
            @test f[9] != ""
        end
    end

    @testset "CSV quotes a comma in the source path" begin
        cf = write_table(joinpath(dir, "a,b.csv"), X)
        r = run_cli("analyze", cf, "--devs", "adev", "--no-ci")
        @test r.code == 0
        @test occursin("\"$cf\"", r.out)
    end

    @testset "--m are averaging factors (τ = m·τ₀)" begin
        r = run_cli("analyze", pf, "--tau0", "2", "--devs", "adev", "--m", "1,2,4", "--no-ci")
        @test r.code == 0
        @test getindex.(csv_rows(r.out), 4) == ["2.0", "4.0", "8.0"]
        @test occursin("\"tau_mode\":\"m\",\"m\":[1,2,4]",
                       run_cli("analyze", pf, "--devs", "adev", "--m", "4,2,1", "--format", "json").out)
    end

    @testset "--grid accepts both spellings" begin
        a = run_cli("analyze", pf, "--devs", "adev", "--grid", "half-octave", "--no-ci").out
        b = run_cli("analyze", pf, "--devs", "adev", "--grid", "halfoctave", "--no-ci").out
        @test a == b
        @test length(csv_rows(a)) > length(csv_rows(run_cli("analyze", pf, "--devs", "adev", "--no-ci").out))
    end

    @testset "reader options" begin
        vf = joinpath(dir, "vals.txt")
        open(vf, "w") do io
            foreach(v -> println(io, v), X)
        end
        r = run_cli("analyze", vf, "--time-col", "0", "--value-col", "1", "--tau0", "1",
                    "--scaling", "1e9", "--detrend", "mean", "--devs", "adev", "--no-ci")
        @test r.code == 0
        sf = joinpath(dir, "semi.txt")
        open(sf, "w") do io
            foreach(((i, v),) -> println(io, i - 1, ";", v), enumerate(X))
        end
        r = run_cli("analyze", sf, "--delim", ";", "--devs", "adev", "--no-ci")
        @test r.code == 0
        @test r.out == replace(run_cli("analyze", pf, "--devs", "adev", "--no-ci").out, pf => sf)
    end

    @testset "reference column" begin
        rf = joinpath(dir, "three.csv")
        open(rf, "w") do io
            for (i, v) in enumerate(X)
                println(io, i - 1, ",", v + 1e-9 * i, ",", 1e-9 * i)
            end
        end
        a = csv_rows(run_cli("analyze", rf, "--reference-col", "3", "--devs", "adev", "--no-ci").out)
        b = csv_rows(run_cli("analyze", pf, "--devs", "adev", "--no-ci").out)
        @test isapprox(parse.(Float64, getindex.(a, 5)), parse.(Float64, getindex.(b, 5)); rtol=1e-6)
        @test run_cli("analyze", rf, "--reference-col", "9").code == 2
        @test run_cli("analyze", rf, "--reference-col", "2").code == 1
    end

    @testset "check" begin
        r = run_cli("check", pf)
        @test r.code == 0
        lines = split(strip(r.out), '\n')
        @test lines[1] == "source,quantity,n,tau0,span_s,gaps,outliers,drift_per_day"
        @test split(lines[2], ',')[3] == "256"
        r = run_cli("check", pf, "--format", "json", "--series", "4000")
        @test r.code == 0
        @test occursin("\"command\":\"check\"", r.out)
        @test occursin("\"n\":256", r.out) && occursin("\"gaps\":0", r.out)
        @test occursin("\"series\":{\"t\":[0.0,1.0,", r.out)
        @test !occursin("\"results\"", r.out)
    end

    @testset "spectrum" begin
        r = run_cli("spectrum", pf)
        @test r.code == 0
        rows = csv_rows(r.out)
        @test split(r.out, '\n')[1] == "f,sy,sx"
        @test length(rows) == 64 ÷ 2       # nperseg 64 = largest power of 2 ≤ 256/4, DC dropped
        r = run_cli("spectrum", pf, "--nperseg", "32", "--format", "json")
        @test occursin("\"spectrum\":{\"nperseg\":32,\"noverlap\":16,", r.out)
        sy = psd_y(PhaseData(X, 1.0); nperseg=32)
        @test parse(Float64, csv_rows(run_cli("spectrum", pf, "--nperseg", "32").out)[1][2]) == sy.psd[2]
    end

    @testset "dynamic" begin
        r = run_cli("dynamic", pf, "--window", "64", "--format", "json")
        @test r.code == 0
        @test occursin("\"dynamic\":{\"deviation\":\"dadev\",\"window\":64,\"step\":32,", r.out)
        ref = dadev(PhaseData(X, 1.0), tau_values(Octave, 64, :adev); window=64, step=32)
        rows = csv_rows(run_cli("dynamic", pf, "--window", "64").out)
        @test length(rows) == length(ref.t) * length(ref.tau)
        @test parse(Float64, rows[1][3]) == ref.dev[1, 1]
        r = run_cli("dynamic", pf, "--dev", "dhdev", "--format", "json")
        @test occursin("\"deviation\":\"dhdev\",\"window\":32,", r.out)
    end

    @testset "synthetic record" begin
        r = run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev", "--format", "json")
        @test r.code == 0
        @test occursin("\"source\":\"synthetic\"", r.out)
        @test run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev").out ==
              run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev").out
    end

    @testset "usage errors -> exit 1" begin
        @test run_cli("analyze").code == 1
        r = run_cli("analyze", pf, "--devs", "bogus")
        @test r.code == 1
        @test occursin("bogus", r.err) && occursin("mhtotdev", r.err)
        r = run_cli("analyze", pf, "--devs", "adev,dadev")
        @test r.code == 1
        @test occursin("sigmatau dynamic", r.err)
        @test run_cli("analyze", pf, "--grid", "weird").code == 1
        @test run_cli("analyze", pf, "--m", "0,2").code == 1
        @test run_cli("analyze", pf, "--phase", "--freq").code == 1
        @test run_cli("analyze", pf, "--confidence", "1.5").code == 1
        @test run_cli("analyze", pf, "--spectrum").code == 1        # CSV has no room for it
        @test run_cli("analyze", pf, pf).code == 1
        @test run_cli("analyze", pf, "--synth", "64").code == 1
        @test run_cli("analyze", pf, "--tau0").code == 1
        @test run_cli("analyze", pf, "--no-ci=yes").code == 1
        @test run_cli("analyze", pf, "--time-col", "-1", "--tau0", "1").code == 1
        @test run_cli("analyze", pf, "--value-col", "0").code == 1
        @test run_cli("analyze", pf, "--tau0", "Inf").code == 1
        @test run_cli("analyze", pf, "--tau0", "-1").code == 1
        @test run_cli("analyze", pf, "--scaling", "NaN").code == 1
    end

    @testset "errors in JSON when asked for JSON" begin
        r = run_cli("analyze", pf, "--devs", "bogus", "--format", "json")
        @test r.code == 1
        @test startswith(r.out, "{\"schema_version\":1,")
        @test occursin("\"error\":\"unknown deviation", r.out)
        @test occursin("\"exit_code\":1}", r.out)
        r = run_cli("analyze", joinpath(dir, "nope.csv"), "--format", "json")
        @test r.code == 2
        @test occursin("\"exit_code\":2", r.out)
    end

    @testset "input error -> exit 2" begin
        r = run_cli("analyze", joinpath(dir, "nope.csv"))
        @test r.code == 2
        @test occursin("no such file", r.err)
        bad = joinpath(dir, "bad.csv")
        write(bad, "0,abc\n1,def\n")
        @test run_cli("analyze", bad).code == 2
    end

    @testset "analysis error -> exit 3" begin
        sf = write_table(joinpath(dir, "short.csv"), X[1:4])
        r = run_cli("analyze", sf, "--tau0", "1", "--devs", "adev")
        @test r.code == 3
        @test occursin("at least 8 samples", r.err)
        gf = joinpath(dir, "gappy.csv")
        open(gf, "w") do io
            for (i, v) in enumerate(X)
                println(io, i - 1, ",", i == 50 ? "NaN" : string(v))
            end
        end
        r = run_cli("spectrum", gf)
        @test r.code == 3
        @test occursin("--fillgaps", r.err)
    end

    @testset "output error -> exit 4" begin
        r = run_cli("analyze", pf, "--devs", "adev", "--out", joinpath(dir, "no", "sub", "out.csv"))
        @test r.code == 4
    end

    @testset "broken pipe on stdout -> clean exit" begin
        e = IOBuffer()
        code = SigmaTauCLI.sigmatau(["analyze", pf, "--devs", "adev"]; out=BrokenPipeIO(), err=e)
        @test code == 0
        @test isempty(String(take!(e)))
    end

    @testset "--out writes what stdout would show" begin
        op = joinpath(dir, "out.csv")
        r = run_cli("analyze", pf, "--devs", "adev", "--out", op)
        @test r.code == 0
        @test occursin("wrote", r.err)
        @test read(op, String) == run_cli("analyze", pf, "--devs", "adev").out
        @test isempty(run_cli("analyze", pf, "--devs", "adev", "--out", op, "--quiet").err)
    end

    exe = get(ENV, "SIGMATAU_EXE", "")
    if !isempty(exe)
        @testset "compiled binary matches" begin
            # elapsed_ms differs between runs; compare everything before it.
            strip_t(s) = replace(s, r"\"elapsed_ms\":[^}]*}" => "")
            # The compiled code may round the last bit differently from the JIT
            # (another CPU target), so the text must match with every number
            # replaced, and the numbers must agree to 1e-12 relative.
            numre = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
            function same_output(a, b)
                replace(a, numre => "#") == replace(b, numre => "#") || return false
                xs = [parse(Float64, m.match) for m in eachmatch(numre, a)]
                ys = [parse(Float64, m.match) for m in eachmatch(numre, b)]
                return all(isapprox(x, y; rtol=1e-12, atol=0) for (x, y) in zip(xs, ys))
            end
            for args in (["analyze", pf, "--devs", "adev,mdev,tdev,totdev,theo1"],
                         ["analyze", pf, "--format", "json", "--spectrum", "--series", "100"],
                         ["check", pf, "--format", "json"],
                         ["spectrum", pf],
                         ["dynamic", pf, "--dev", "dhdev"],
                         ["analyze", pf, "--devs", "bogus", "--format", "json"],
                         ["completion", "bash"],
                         ["version"])
                o = IOBuffer()
                p = run(pipeline(ignorestatus(`$exe $args`); stdout=o, stderr=devnull))
                want = run_cli(args...)
                @test same_output(strip_t(String(take!(o))), strip_t(want.out))
                @test p.exitcode == want.code
            end
        end
    end
end

end
