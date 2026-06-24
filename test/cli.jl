# test/cli.jl — non-interactive CLI (`SigmaTau.cli_main`) and `bin/sigmatau.jl`.
#
# Every in-process case calls `cli_main` with injected IOBuffers so exit codes
# and stdout/stderr can be asserted without touching the real streams. Data is
# deterministic (no Random) so results are reproducible.

using Test
using SigmaTau

# Run cli_main with captured streams. Returns (code, out, err).
function run_cli(args...)
    o = IOBuffer()
    e = IOBuffer()
    code = SigmaTau.cli_main(collect(String, args); out = o, err = e)
    return (code = code, out = String(take!(o)), err = String(take!(e)))
end

# Deterministic "phase" record (seconds), long enough for octave grids.
const _CLI_X = cumsum(sinpi.((1:128) ./ 7.3)) ./ 1e9

# An output stream that always fails with an EPIPE IOError, simulating a
# downstream consumer (e.g. `| head`) that closed the pipe.
struct BrokenPipeIO <: IO end
_epipe() = Base.IOError("write: broken pipe (EPIPE)", Base.UV_EPIPE)
Base.write(::BrokenPipeIO, ::UInt8) = throw(_epipe())
Base.unsafe_write(::BrokenPipeIO, ::Ptr{UInt8}, ::UInt) = throw(_epipe())

# Write a 2-column time,value table (extension picks the delimiter for the reader).
function _write_table(path, x; tau0 = 1.0)
    open(path, "w") do io
        for (i, v) in enumerate(x)
            println(io, (i - 1) * tau0, ",", v)
        end
    end
    return path
end

@testset "version" begin
    for v in (["version"], ["--version"], ["-V"])
        r = run_cli(v...)
        @test r.code == 0
        @test occursin(string(pkgversion(SigmaTau)), r.out)
    end
end

@testset "help" begin
    for h in (["help"], ["--help"], ["-h"], String[])
        r = run_cli(h...)
        @test r.code == 0
        @test occursin("analyze", r.out)
        @test occursin("completion", r.out)
        @test occursin("--devs", r.out)
        @test occursin("--out", r.out)
    end
end

@testset "unknown command / option" begin
    r = run_cli("frobnicate")
    @test r.code == 1
    @test occursin("unknown command", r.err)

    r2 = run_cli("--bogus")
    @test r2.code == 1
    @test occursin("unknown option", r2.err)
end

@testset "completion scripts" begin
    for sh in ("bash", "zsh", "fish", "powershell")
        r = run_cli("completion", sh)
        @test r.code == 0
        s = r.out
        # fish names options as `-l devs`; the others as `--devs`.
        optpat = sh == "fish" ? "-l " : "--"
        @test occursin("analyze", s)
        @test occursin("completion", s)
        @test occursin(optpat * "devs", s)
        @test occursin(optpat * "out", s)
        @test occursin("dadev", s) && occursin("dhdev", s)   # dynamic devs in vocabulary
        @test occursin("octave", s) && occursin("decade", s)
        @test occursin("csv", s) && occursin("json", s)
        @test occursin("true", s) && occursin("false", s)
    end
    # Per-shell native file-completion directive must be present.
    @test occursin("-o default", run_cli("completion", "bash").out)
    @test occursin("_files", run_cli("completion", "zsh").out)
    fishout = run_cli("completion", "fish").out
    @test occursin("-F", fishout)                           # force-files for --out
    @test occursin("'__fish_seen_subcommand_from analyze' -l devs", fishout)  # value rule w/o -f
    @test occursin("Register-ArgumentCompleter", run_cli("completion", "powershell").out)

    # Unknown / missing shell.
    r = run_cli("completion", "tcsh")
    @test r.code == 1
    @test occursin("tcsh", r.err)
    @test run_cli("completion").code == 1
end

@testset "analyze --help" begin
    r = run_cli("analyze", "--help")
    @test r.code == 0
    @test occursin("analyze FILE", r.out)
    @test occursin("--devs", r.out)
end

mktempdir() do dir
    pf = _write_table(joinpath(dir, "phase.csv"), _CLI_X)

    @testset "analyze phase -> CSV (stdout)" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev,mdev",
                    "--tau", "octave", "--format", "csv")
        @test r.code == 0
        lines = split(strip(r.out), '\n')
        @test lines[1] == "source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type"
        @test length(lines) > 1
        @test all(l -> length(split(l, ',')) == 10, lines[2:end])
        @test any(l -> occursin(",adev,", l), lines[2:end])
        @test any(l -> occursin(",mdev,", l), lines[2:end])
        @test all(l -> occursin(",phase,", l), lines[2:end])
    end

    @testset "analyze frequency -> CSV" begin
        ff = _write_table(joinpath(dir, "freq.csv"), diff(_CLI_X))
        r = run_cli("analyze", ff, "--freq", "--tau0", "1", "--devs", "adev", "--format", "csv")
        @test r.code == 0
        @test all(l -> occursin(",frequency,", l), split(strip(r.out), '\n')[2:end])
    end

    @testset "analyze -> JSON" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev", "--format", "json")
        @test r.code == 0
        for k in ("\"schema_version\"", "\"analysis\"", "\"results\"", "\"points\"",
                  "\"deviations\"", "\"noise_type\"", "\"tau_mode\"")
            @test occursin(k, r.out)
        end
    end

    @testset "JSON ci=false -> nulls, never NaN" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev",
                    "--format", "json", "--ci", "false")
        @test r.code == 0
        @test occursin("\"confidence\":", r.out)
        @test occursin("null", r.out)
        @test !occursin("NaN", r.out)
    end

    @testset "CSV ci=false blanks CI columns, never NaN" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev",
                    "--format", "csv", "--ci", "false")
        @test r.code == 0
        @test !occursin("NaN", r.out)
        for l in split(strip(r.out), '\n')[2:end]
            f = split(l, ',')
            @test f[6] == "" && f[7] == "" && f[8] == "" && f[10] == ""  # ci_lower/upper/edf/noise
            @test f[9] != ""                                             # neff always present
        end
    end

    @testset "CSV quotes a comma in the source path" begin
        cf = _write_table(joinpath(dir, "a,b.csv"), _CLI_X)
        r = run_cli("analyze", cf, "--tau0", "1", "--devs", "adev", "--ci", "false")
        @test r.code == 0
        lines = split(strip(r.out), '\n')
        @test startswith(lines[2], "\"")
        @test occursin("\"$cf\"", r.out)
    end

    @testset "explicit --taus are averaging factors (tau = m*tau0)" begin
        r = run_cli("analyze", pf, "--tau0", "2", "--devs", "adev",
                    "--taus", "1,2,4", "--ci", "false", "--format", "csv")
        @test r.code == 0
        taus_out = [split(l, ',')[4] for l in split(strip(r.out), '\n')[2:end]]
        @test taus_out == ["2.0", "4.0", "8.0"]
    end

    @testset "reader pass-through kwargs" begin
        vf = joinpath(dir, "vals.txt")
        open(vf, "w") do io
            for v in _CLI_X
                println(io, v)
            end
        end
        r = run_cli("analyze", vf, "--time-col", "0", "--value-col", "1", "--tau0", "1",
                    "--scaling", "1e9", "--detrend", "mean", "--devs", "adev", "--ci", "false")
        @test r.code == 0
    end

    @testset "usage errors -> exit 1" begin
        @test run_cli("analyze").code == 1                                   # missing FILE
        @test run_cli("analyze", pf, "--devs", "bogus").code == 1            # unknown dev
        let r = run_cli("analyze", pf, "--devs", "bogus")
            @test occursin("bogus", r.err)
            @test occursin("adev", r.err)                                    # lists valid
            @test !occursin("MethodError", r.err)                           # no raw trace
        end
        let r = run_cli("analyze", pf, "--devs", "adev,dadev")               # dynamic dev
            @test r.code == 1
            @test occursin("dadev", r.err) && occursin("dynamic", r.err)
        end
        @test run_cli("analyze", pf, "--tau", "octave", "--taus", "1,2").code == 1
        @test run_cli("analyze", pf, "--tau", "weird").code == 1
        @test run_cli("analyze", pf, "--ci", "maybe").code == 1
        @test run_cli("analyze", pf, "--phase", "--freq").code == 1
    end

    @testset "input error -> exit 2" begin
        r = run_cli("analyze", joinpath(dir, "nope.csv"), "--tau0", "1")
        @test r.code == 2
        @test !occursin("MethodError", r.err)
    end

    @testset "analysis error -> exit 3" begin
        sf = _write_table(joinpath(dir, "short.csv"), _CLI_X[1:2])
        r = run_cli("analyze", sf, "--tau0", "1", "--devs", "adev")
        @test r.code == 3
    end

    @testset "output error -> exit 4" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev",
                    "--out", joinpath(dir, "no", "sub", "out.csv"))
        @test r.code == 4
    end

    @testset "broken pipe on stdout -> clean exit, no error" begin
        e = IOBuffer()
        code = SigmaTau.cli_main(["analyze", pf, "--tau0", "1", "--devs", "adev"];
                                 out = BrokenPipeIO(), err = e)
        @test code == 0
        @test isempty(String(take!(e)))
        @test SigmaTau._is_broken_pipe(_epipe())
        @test !SigmaTau._is_broken_pipe(ErrorException("nope"))
    end

    @testset "valid --out writes a file matching stdout" begin
        op = joinpath(dir, "out.csv")
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev", "--out", op)
        @test r.code == 0
        @test isfile(op)
        content = read(op, String)
        @test startswith(content, "source,quantity,deviation")
        r2 = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev")
        @test strip(content) == strip(r2.out)
    end
end

@testset "bin/sigmatau.jl version (subprocess)" begin
    proj = dirname(@__DIR__)                       # test/ -> repo root
    binscript = joinpath(proj, "bin", "sigmatau.jl")
    @test isfile(binscript)
    cmd = `$(Base.julia_cmd()) --project=$proj $binscript version`
    out = read(ignorestatus(cmd), String)
    @test occursin(string(pkgversion(SigmaTau)), out)
end
