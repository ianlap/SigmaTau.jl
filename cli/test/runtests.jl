# Tests for the sigmatau command line, run in process:
#
#   julia --project=cli cli/test/runtests.jl
#
# With SIGMATAU_EXE set to a compiled binary (cli/build/bin/sigmatau), the
# last testset also runs it and checks that it answers exactly like the
# in-process call. The completion scripts run in each of bash, zsh, fish and
# pwsh that is installed; CI has all four.

using Test
using SigmaTau
using SigmaTauCLI
using Random

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

# Each shell loads its completion script and prints, one line per command
# line (current word last), the words it offers. zsh completes only at a
# terminal, so its completion builtins are stubbed.
quoted(w) = join(repr.(w), ' ')
const SHELL_DRIVERS = Dict(
    "bash" => (`bash`, (s, q) -> "source $s\n" * join(("COMP_WORDS=($(quoted(w))); COMP_CWORD=$(length(w) - 1); COMPREPLY=(); _sigmatau; echo \"\${COMPREPLY[*]}\"" for w in q), "\n")),
    "zsh" => (`zsh`, (s, q) -> """
        _arguments() { if (( CURRENT == 2 )); then state=cmd; else state=args; shift words; (( CURRENT-- )); fi }
        _values() { [[ \$1 == -s ]] && shift 2; shift; print -rn -- "\$* " }
        compadd() { print -rn -- "\${(P)2} " }
        _files() { print -rn -- "FILES " }
        words=(sigmatau ''); CURRENT=2; source $s >/dev/null
        """ * join(("words=($(quoted(w))); CURRENT=$(length(w)); PREFIX=$(repr(w[end])); _sigmatau; print" for w in q), "\n")),
    "fish" => (`fish`, (s, q) -> "source $s\n" * join(("echo (complete -C $(repr(join(w, ' '))) | string replace -r '\\t.*' '')" for w in q), "\n")),
    "powershell" => (`pwsh -NoProfile`, (s, q) -> ". $s\n" * join(("\$q = '$(join(w, ' '))'; (TabExpansion2 -inputScript \$q -cursorColumn \$q.Length).CompletionMatches.CompletionText -join ' '" for w in q), "\n")),
)

function shell_completions(sh, queries)
    exe, driver = SHELL_DRIVERS[sh]
    mktempdir() do dir
        ext = sh == "powershell" ? ".ps1" : ""
        script, drv = joinpath(dir, "sigmatau" * ext), joinpath(dir, "driver" * ext)
        write(script, run_cli("completion", sh).out)
        write(drv, driver(script, queries))
        # A `sigmatau` on the PATH for PowerShell to find, run in a folder with no files to offer.
        stub = joinpath(mkdir(joinpath(dir, "bin")), "sigmatau")
        write(stub, "#!/bin/sh\n"); chmod(stub, 0o755)
        cmd = addenv(`$exe $drv`, "PATH" => dirname(stub) * ":" * ENV["PATH"])
        out = cd(() -> read(cmd, String), mkdir(joinpath(dir, "empty")))
        return [Set(split(l)) for l in split(chomp(out), '\n')]
    end
end

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
        for c in ("analyze", "check", "spectrum", "dynamic", "live", "completion")
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
    # A signature is set apart from its description, and σ keeps its case.
    @test all(d -> occursin("  " * d * "\n", run_cli("help").out), last.(SigmaTauCLI.COMMAND_HELP))
    @test occursin("\nThe σ(τ)", run_cli("help", "model").out)
    @test !occursin("FILE", run_cli("help", "live").out)
end

@testset "the reference names every command and option" begin
    page = read(joinpath(@__DIR__, "..", "..", "docs", "src", "reference", "cli.md"), String)
    for (c, arg, _) in SigmaTauCLI.COMMAND_HELP
        @test occursin("| `" * (isempty(arg) ? c : c * " " * arg) * "` |", page)
    end
    options = Set(op.long for op in SigmaTauCLI.OPTS)
    named = Set(m.match for m in eachmatch(r"--[a-z][a-z0-9-]*", page))
    @test setdiff(options, named) == Set{String}()      # every option is on the page
    @test setdiff(named, options) == Set{String}()      # and the page names no other
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
    # Every command offers its options, and every option with fixed values offers them.
    queries = [["sigmatau", ""], ["sigmatau", "completion", ""]]
    want = [Set(SigmaTauCLI.subcommands()), Set(SigmaTauCLI.SHELLS)]
    for c in SigmaTauCLI.COMMANDS
        push!(queries, ["sigmatau", c, "--"]); push!(want, Set(SigmaTauCLI.cmd_options(c)))
    end
    for op in SigmaTauCLI.valued_options()
        push!(queries, ["sigmatau", op.cmds[1], op.long, ""]); push!(want, Set(SigmaTauCLI.completion_values(op)))
    end
    shells = filter(sh -> Sys.which(SHELL_DRIVERS[sh][1].exec[1]) !== nothing, SigmaTauCLI.SHELLS)
    haskey(ENV, "GITHUB_ACTIONS") && @test length(shells) == length(SigmaTauCLI.SHELLS)
    for sh in shells
        @test shell_completions(sh, queries) == want
    end
    # A comma list completes its last item.
    for sh in intersect(shells, ["bash", "fish"])
        @test shell_completions(sh, [["sigmatau", "analyze", "--devs", "adev,mht"]]) == [Set(["adev,mhtotdev"])]
    end
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
        for mode in (:published, :stable32)
            ref = htotdev(pd, [3, 9]; correct_bias=mode)
            rows = csv_rows(run_cli("analyze", pf, "--devs", "htotdev", "--m", "3,9", "--bias", String(mode)).out)
            @test parse.(Float64, getindex.(rows, 5)) == ref.dev
        end
        @test run_cli("analyze", pf, "--bias", "sp1065").code == 1
        for (word, noise) in (("FLFM", :FLFM), ("fit", powerlaw_fit(pd)))
            ref = mtotdev(pd, [3, 9]; noise=noise)
            rows = csv_rows(run_cli("analyze", pf, "--devs", "mtotdev", "--m", "3,9", "--noise", word).out)
            @test parse.(Float64, getindex.(rows, 5)) == ref.dev
            @test parse.(Float64, getindex.(rows, 8)) == ref.edf
        end
        @test run_cli("analyze", pf, "--noise", "pink").code == 1
    end

    @testset "frequency input" begin
        ff = write_table(joinpath(dir, "freq.csv"), diff(X))
        r = run_cli("analyze", ff, "--freq", "--tau0", "1", "--devs", "adev")
        @test r.code == 0
        @test all(f -> f[2] == "frequency", csv_rows(r.out))
        # Every frequency sample is used, as in the library's FrequencyData methods.
        r = run_cli("analyze", ff, "--freq", "--tau0", "1", "--devs", "adev", "--m", "1,2,4")
        ref = adev(FrequencyData(diff(X), 1.0), [1, 2, 4]; ci=false)
        @test parse.(Float64, getindex.(csv_rows(r.out), 5)) ≈ ref.dev rtol=1e-12
    end

    @testset "analyze -> JSON" begin
        r = run_cli("analyze", pf, "--tau0", "1", "--devs", "adev,tdev", "--format", "json")
        @test r.code == 0
        for k in ("\"schema_version\":1", "\"command\":\"analyze\"", "\"analysis\"", "\"results\"",
                  "\"points\"", "\"deviations\"", "\"noise_type\"", "\"tau_mode\":\"octave\"",
                  "\"axis\":\"sigma_x\"", "\"label\":\"TDEV\"", "\"record\"", "\"elapsed_ms\"",
                  "\"bias_correction\":\"exact\"", "\"noise\":\"auto\"")
            @test occursin(k, r.out)
        end
        @test !occursin("\"spectrum\"", r.out)
        @test !occursin("\"series\"", r.out)
        r = run_cli("analyze", pf, "--format=json", "--spectrum", "--series", "10")
        @test r.code == 0
        @test occursin("\"spectrum\":{\"nperseg\":64", r.out)
        @test occursin("\"series\":{\"t\":[", r.out)
        # thinned by stretch: an outlier and a gap between strided samples stay
        x = zeros(1000); x[437] = 1.0; x[800] = NaN
        idx = SigmaTauCLI.series_index(x, 30)
        @test 437 in idx && 800 in idx && issorted(idx) && length(idx) <= 30
        @test idx[1] == 1 && idx[end] == 1000
        @test SigmaTauCLI.series_index(x[1:20], 30) == 1:20
        r = run_cli("analyze", pf, "--format=json", "--devs", "adev,tdev,totdev", "--fit")
        @test r.code == 0
        @test occursin("\"powerlaw\":{\"deviation\":\"mhdev\",\"alphas\":[2,1,0,-1,-2]", r.out)
        @test occursin("\"curves\":[{\"deviation\":\"adev\"", r.out)
        @test occursin("{\"deviation\":\"tdev\"", r.out) && !occursin("{\"deviation\":\"totdev\",\"tau\"", r.out)
        @test occursin("\"sigma_at\":{\"tau\":[1.0,10.0,100.0]", r.out)
        @test run_cli("analyze", pf, "--fit").code != 0
        r = run_cli("analyze", "--synth", "4096", "--format=json", "--devs", "mhdev", "--component", "WHFM::64", "--component", "-2:256:")
        @test r.code == 0
        @test occursin("\"alphas\":[0,-2],\"h\":[", r.out) && occursin("\"mode\":\"components\"", r.out)
        @test occursin("\"components\":[{\"alpha\":0,", r.out) && occursin("\"residual\":{\"tau\":[", r.out)
        @test occursin("\"clock_model\":{\"R\":0.0,\"q1\":", r.out) && occursin("\"q2\":", r.out)
        @test occursin("\"error\":", run_cli("analyze", pf, "--format=json", "--component", "WHFM:1e9:").out)
        @test run_cli("analyze", pf, "--format=json", "--component", "PINK:1:2").code == 1
        r = run_cli("analyze", pf, "--format=json", "--devs", "adev,tdev,mtie", "--mask", "G8262_opt1")
        @test r.code == 0 && occursin("\"masks\":[{\"name\":\"G8262_opt1\"", r.out)
        @test occursin("\"checks\":[{\"deviation\":\"tdev\"", r.out) && occursin("{\"deviation\":\"mtie\"", r.out)
        @test run_cli("analyze", pf, "--mask", "G8262_opt1").code != 0      # JSON only
        @test run_cli("analyze", pf, "--format=json", "--mask", "G999").code != 0
        @test count("\"title\":", run_cli("analyze", pf, "--format=json", "--mask", "all").out) == 9
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

    @testset "instrument formats" begin
        fx(f) = joinpath(@__DIR__, "..", "..", "test", "fixtures", "formats", f)
        tim = fx("timelab.tim")
        a = csv_rows(run_cli("analyze", tim, "--devs", "adev", "--no-ci").out)
        b = adev(read_data(tim); ci=false)
        @test parse.(Float64, getindex.(a, 5)) ≈ b.dev
        @test parse(Float64, a[1][4]) == 0.1                     # τ₀ from the header
        @test occursin("frequency", run_cli("analyze", fx("counter_53230a_freq.csv"), "--tau0", "1",
                                            "--devs", "adev", "--no-ci", "--fillgaps").out)
        # A missing frequency reading (9.91E+37) must be filled, not integrated as NaN.
        @test run_cli("analyze", fx("counter_53230a_freq.csv"), "--tau0", "1").code == 2
        @test occursin("\"start\":\"2026-11-08T12:00:00Z\",\"gaps\":1",
                       run_cli("check", fx("gpsdo.csv"), "--format", "json").out)
        @test run_cli("analyze", fx("counter_53230a_freq.csv")).code == 2   # no gate time
        g = run_cli("check", fx("gpsdo.csv"), "--column", "PPS offset (ns)", "--format", "json")
        @test g.code == 0 && occursin("\"quantity\":\"phase\"", g.out)
        o = run_cli("check", fx("gpsdo.csv"), "--column", "Temperature (C)", "--phase", "--format", "json")
        @test o.code == 0 && occursin("\"quantity\":\"phase\"", o.out)
        @test run_cli("check", fx("gpsdo.csv"), "--column", "OSC offset (ppb)", "--phase").code == 2
        @test run_cli("analyze", tim, "--input-format", "counter", "--tau0", "1").code == 2
        # A BIPM UTC - UTC(k) file is recognized by its header: phase, five days apart, two rows missing.
        b = run_cli("check", fx("bipm_utc_lab.txt"), "--format", "json")
        @test b.code == 0 && occursin("\"quantity\":\"phase\"", b.out) && occursin("\"tau0\":432000.0", b.out)
        @test occursin("\"gaps\":2", b.out)
        # A RINEX clock file: the clock named by --column, phase every 30 s.
        r = run_cli("check", fx("rinex_clk_synth.clk"), "--column", "E02", "--format", "json")
        @test r.code == 0 && occursin("\"quantity\":\"phase\"", r.out) && occursin("\"tau0\":30.0", r.out)
        @test run_cli("check", fx("rinex_clk_synth.clk"), "--column", "G99").code == 2
        @test run_cli("analyze", tim, "--reference-col", "2").code == 1
        semi = joinpath(dir, "semi.txt")
        write(semi, "time;phase (ns)\n" * join(string.(0:63, ";", X[1:64] .* 1e9), "\n"))
        d = csv_rows(run_cli("analyze", semi, "--input-format", "table", "--delim", ";",
                             "--devs", "adev", "--no-ci").out)
        @test parse.(Float64, getindex.(d, 5)) ≈ adev(PhaseData(X[1:64], 1.0); ci=false).dev rtol=1e-9
        # The filled second of the GPSDO log stays out of the drift fit.
        o = SigmaTauCLI.Options()
        o.input, o.fillgaps = fx("gpsdo.csv"), true
        SigmaTauCLI.read_record(o)
        @test findall(!, o.measured) == [24]
    end

    @testset "cross" begin
        X2 = X .+ 1e-11 .* sinpi.((1:256) ./ 3.1)
        f2 = write_table(joinpath(dir, "second.csv"), X2)
        r = run_cli("cross", pf, "--with", f2, "--m", "1,2,4")
        @test r.code == 0
        rows = csv_rows(r.out)
        lib = cross_adev(PhaseData(X, 1.0), PhaseData(X2, 1.0), [1, 2, 4])
        @test parse.(Float64, getindex.(rows, 2)) ≈ lib.covar
        @test occursin("\"cross_spectrum\"", run_cli("cross", pf, "--with", f2, "--format", "json").out)
        # --remove-drift reports the first channel's fit, not the second's.
        j = run_cli("cross", pf, "--with", f2, "--remove-drift", "linear", "--format", "json").out
        c = run_cli("check", pf, "--remove-drift", "linear", "--format", "json").out
        offset(s) = match(r"\"removed_offset\":([^,}]+)", s)[1]
        @test offset(j) == offset(c)
        @test run_cli("cross", pf).code == 1
        @test run_cli("cross", pf, "--with", joinpath(dir, "missing.csv")).code == 2
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

    @testset "timescale" begin
        # Three clocks in columns 2–4 against one reference, the third twice as noisy.
        cc, _ = simulate_ensemble([ClockModel(; tau0=1.0, q1=q) for q in (1e-24, 1e-24, 4e-24)], 4096;
                                  reference=0, rng=Random.Xoshiro(3))
        tf = joinpath(dir, "clocks.csv")
        open(tf, "w") do io
            for k in 1:4096
                println(io, k - 1, ",", join(cc.x[k, :], ","))
            end
        end
        rows = csv_rows(run_cli("timescale", tf, "--value-cols", "2,3,4").out)
        @test getindex.(rows, 1) == ["column 2", "column 3", "column 4"]
        w = parse.(Float64, getindex.(rows, 3))
        @test sum(w) ≈ 1 && w[3] < w[1] / 2
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "ewfa", "--format=json", "--series", "8")
        @test r.code == 0
        @test occursin("\"algorithm\":\"ewfa\"", r.out) && occursin("\"weight\":[0.3333333333333333,", r.out)
        @test occursin("\"adev\":{\"tau\":[1.0,", r.out) && occursin("\"series\":{\"t\":[0.0,", r.out)
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "at2", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"algorithm\":\"at2\"", r.out)
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "at2", "--steps", "real", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"steps\":\"real\"", r.out)
        @test run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "at1", "--steps", "best").code != 0
        # The Kalman scales fit a model per clock; --pivot picks the clock they read at.
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "greenhall", "--pivot", "2", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"algorithm\":\"greenhall\"", r.out) && occursin("\"pivot\":2", r.out)
        @test run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "jones_tryon", "--format=json", "--series", "0").code == 0
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "stein", "--wmax", "0.5", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"algorithm\":\"stein\"", r.out)
        @test run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "greenhall", "--pivot", "4").code != 0
        # --flicker fits flicker states per clock; these clocks have no flicker, so none.
        r = run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "greenhall", "--flicker", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"flicker\":{\"fm\":[0,0,0],\"pm\":[0,0,0]}", r.out)
        @test run_cli("timescale", tf, "--value-cols", "2,3,4", "--algorithm", "stein", "--flicker").code != 0
        @test run_cli("timescale", tf, "--value-cols", "2,3", "--algorithm", "greenhall").code != 0
        @test run_cli("timescale", "--synth", "2048", "--format=json").code == 0
        # --clocks sets how many clocks are simulated (four by default).
        r = run_cli("timescale", "--synth", "2048", "--clocks", "6", "--format=json", "--series", "0")
        @test r.code == 0 && occursin("\"clocks\":[\"clock 1\",\"clock 2\",\"clock 3\",\"clock 4\",\"clock 5\",\"clock 6\"]", r.out)
        @test length(csv_rows(run_cli("timescale", "--synth", "2048").out)) == 4
        @test run_cli("timescale", tf).code == 1                       # no --value-cols
        @test run_cli("timescale", tf, "--value-cols", "2,2").code == 1
        @test run_cli("timescale", tf, "--value-cols", "2,9").code == 2
        # A table with a line of names is read clock by clock, in the units
        # and on the MJD axis its names give; --header reads the same axis.
        nf = joinpath(dir, "named.csv")
        open(nf, "w") do io
            println(io, "MJD,a (ns),b (ns),c (ns)")
            for k in 1:64
                println(io, 60000 + 5(k - 1), ",", join(cc.x[k, :] .* 1e9, ","))
            end
        end
        for extra in ((), ("--header", "1", "--scaling", "1e-9"))
            r = run_cli("timescale", nf, "--value-cols", "2,3,4", extra..., "--format=json", "--series", "0")
            @test r.code == 0 && occursin("\"tau0\":432000.0,\"n\":64,", r.out)
        end
        @test run_cli("timescale", nf, "--value-cols", "2,3", "--column", "a (ns)").code == 1
        # Errors name the command line's flags, not read_data's keyword.
        r = run_cli("analyze", nf)
        @test r.code == 2 && occursin("pass --column with one of MJD, a (ns)", r.err) && !occursin("column=", r.err)
        r = run_cli("timescale", nf, "--value-cols", "2,7")
        @test r.code == 2 && occursin("column 7 exceeds 4 columns", r.err)
    end

    @testset "drift" begin
        # Frequency with an injected aging of 2e-15 per second.
        y = 1e-11 .* sinpi.((1:400) ./ 13) .+ 2e-15 .* (0:399) .* 60
        ff = write_table(joinpath(dir, "aging.csv"), y; tau0=60.0)
        ref = drift_fit(FrequencyData(y, 60.0))
        rows = csv_rows(run_cli("drift", ff, "--freq").out)
        @test rows[1][3] == "linear" && isempty(rows[1][8])
        @test parse(Float64, rows[1][10]) ≈ ref.rate rtol=1e-9
        r = run_cli("drift", ff, "--freq", "--model", "log", "--format", "json", "--series", "50")
        @test r.code == 0 && occursin("\"model\":\"log\"", r.out) && occursin("\"fit\":[", r.out)
        @test run_cli("drift", ff, "--model", "cubic").code == 1
        # Four phase samples: three frequency steps fit a line but not a log.
        sf = write_table(joinpath(dir, "short.csv"), [0.0, 1.0, 3.0, 6.0] .* 1e-9)
        @test run_cli("drift", sf).code == 0
        @test run_cli("drift", sf, "--model", "log").code == 3
    end

    @testset "holdover" begin
        p, _ = simulate_clock(ClockModel(; tau0=1.0, q1=1e-22, q2=1e-27), 4096; rng=Random.Xoshiro(5))
        hf = write_table(joinpath(dir, "free.csv"), p.x)
        q = holdover(p; train=512, tau=[16, 64, 256]).quantile
        ref = holdover(p; train=512, tau=[16, 64, 256], limit=sqrt(q[2] * q[3]))
        rows = csv_rows(run_cli("holdover", hf, "--train", "512", "--tau", "256,16,64", "--limit", string(ref.limit)).out)
        @test parse.(Float64, getindex.(rows, 3)) == ref.tau && parse.(Float64, getindex.(rows, 4)) ≈ ref.rms
        @test parse(Float64, rows[1][10]) ≈ ref.time_to_limit
        @test isempty(csv_rows(run_cli("holdover", hf).out)[1][9])                 # no limit: blank
        r = run_cli("holdover", hf, "--train", "512", "--fit-drift", "--no-anchor", "--limit", "1", "--fit", "--format", "json")
        @test r.code == 0 && occursin("\"drift\":true,\"anchor\":false,\"confidence\":0.95", r.out)
        @test occursin("\"limit\":1.0,\"time_to_limit\":null,\"predicted\":{\"alphas\":[2,1,0,-1,-2]", r.out)
        # A stretch twice each holdover time, beside √6 TDEV at the same τ.
        g = holdover(p; grow=2, tau=[1, 16, 256])
        rows = csv_rows(run_cli("holdover", hf, "--grow", "2", "--tau", "1,16,256").out)
        @test parse.(Float64, getindex.(rows, 2)) == 2 .* g.tau && parse.(Float64, getindex.(rows, 4)) ≈ g.rms
        r = run_cli("holdover", hf, "--grow", "2", "--tau", "1,16,256", "--fit", "--format", "json")
        @test occursin("\"train\":null,\"grow\":2,", r.out) && occursin("\"predicted\":{\"alphas\"", r.out)
        @test parse.(Float64, split(match(r"\"reference\":\[([^\]]*)\]", r.out)[1], ',')) ≈ sqrt(6) .* tdev(p, [1, 16, 256]; ci=false).dev
        @test run_cli("holdover", hf, "--grow", "2", "--train", "512").code == 1
        @test run_cli("holdover", hf, "--grow", "2", "--fit-drift").code == 1
        # A Kalman filter as the predictor, its model learned on the first stretch alone, flicker kept.
        kmodel(x) = clock_model(powerlaw_fit(PhaseData(x[1:512], 1.0)); flicker=true, tau_max=length(x))
        k = holdover(p, kmodel(p.x); train=512, tau=[16, 64, 256])
        rows = csv_rows(run_cli("holdover", hf, "--train", "512", "--tau", "16,64,256", "--kalman").out)
        @test parse.(Float64, getindex.(rows, 4)) ≈ k.rms && parse.(Float64, getindex.(rows, 11)) ≈ k.claimed
        r = run_cli("holdover", hf, "--train", "512", "--kalman", "--format", "json")
        @test r.code == 0 && occursin("\"claimed\":[", r.out) && occursin("\"model\":{\"R\":", r.out)
        ff = noise_gen(PhaseData, 4096, 1.0; h=Dict(-1 => 3e-25), rng=Random.Xoshiro(7))
        r = run_cli("holdover", write_table(joinpath(dir, "flicker.csv"), ff.x), "--train", "512", "--tau", "1,16", "--kalman", "--format", "json")
        @test parse.(Float64, split(match(r"\"claimed\":\[([^\]]*)\]", r.out)[1], ',')) ≈ holdover(ff, kmodel(ff.x); train=512, tau=[1, 16]).claimed
        @test parse(Float64, match(r"\"flicker_fm\":([^,}]*)", r.out)[1]) > 1e-26
        @test run_cli("holdover", hf, "--kalman", "--no-anchor").code == 1
        @test run_cli("holdover", hf, "--kalman", "--grow", "2").code == 1
        @test run_cli("holdover", hf, "--train", "0.5").code == 3
        @test run_cli("holdover", hf, "--limit", "-1").code == 1
        @test run_cli("holdover", hf, "--fit").code == 1                           # CSV has no room for it
    end

    @testset "steer" begin
        p, _ = simulate_clock(ClockModel(; tau0=1.0, R=1e-20, q1=1e-22, q2=1e-26), 8192; rng=Random.Xoshiro(5))
        sf = write_table(joinpath(dir, "steer.csv"), p.x)
        # The model comes from the first stretch alone; every result from its end on.
        model = clock_model(powerlaw_fit(PhaseData(p.x[1:2048], 1.0)); flicker=true, tau_max=8192)
        g = steering_gains(30.0, 1.0)
        s = steer(p, model, g)
        rows = csv_rows(run_cli("steer", sf, "--time-constant", "30", "--m", "1,16,256").out)
        @test parse.(Float64, getindex.(rows, 3)) ≈ adev(PhaseData(p.x[2049:end], 1.0), [1, 16, 256]; ci=false).dev
        @test parse.(Float64, getindex.(rows, 4)) ≈ adev(PhaseData(s.steered.x[2049:end], 1.0), [1, 16, 256]; ci=false).dev
        @test parse.(Float64, getindex.(rows, 5)) ≈ steer(model, g; tau=[1.0, 16.0, 256.0], readings=true).adev
        # Gains from a cost, a steer every fourth sample: no prediction at that interval.
        r = run_cli("steer", sf, "--cost", "1:0:1e4", "--interval", "4", "--train", "1024", "--format", "json", "--series", "50")
        @test r.code == 0 && occursin("\"predicted\":{\"error\":", r.out) && occursin("\"correction\":[", r.out)
        @test parse(Float64, match(r"\"gains\":\{\"phase\":([^,]*)", r.out)[1]) ≈ lqg_gains(4.0; phase=1, steer=1e4).phase
        @test all(isempty, getindex.(csv_rows(run_cli("steer", sf, "--cost", "1:0:1e4", "--interval", "4").out), 5))
        r = run_cli("steer", sf, "--time-constant", "30", "--integral", "--fit-drift", "--format", "json")
        @test r.code == 0 && occursin("\"drift\":true", r.out) && !occursin("\"integral\":0,", r.out)
        @test run_cli("steer", sf).code == 1
        @test run_cli("steer", sf, "--time-constant", "30", "--cost", "1:0:1").code == 1
        @test run_cli("steer", sf, "--cost", "1:0:1", "--integral").code == 1
        @test run_cli("steer", sf, "--cost", "0:0:1").code == 1
        @test run_cli("steer", sf, "--time-constant", "30", "--interval", "1.5").code == 3
        @test run_cli("steer", sf, "--time-constant", "30", "--train", "8190").code == 3
    end

    @testset "--remove-drift" begin
        # A known offset and drift on top of X: linear removes the offset,
        # quadratic both, and what was removed is reported.
        t = 0:length(X)-1
        df = write_table(joinpath(dir, "drift.csv"), X .+ 1e-9 .+ 3e-11 .* t .+ 2e-12 .* t .^ 2)
        r = run_cli("check", df, "--remove-drift", "quadratic", "--format", "json")
        @test r.code == 0
        @test occursin("\"remove_drift\":\"quadratic\"", r.out)
        off = parse(Float64, match(r"\"removed_offset\":([^,}]+)", r.out)[1])
        dr = parse(Float64, match(r"\"removed_drift_per_day\":([^,}]+)", r.out)[1])
        @test isapprox(off, 3e-11 + 2 * 2e-12 * (length(X) - 1) / 2; rtol=0.02)
        @test isapprox(dr, 2 * 2e-12 * 86400; rtol=0.02)
        # Removing a parabola leaves HDEV, which ignores drift, nearly alone
        # and matches the library's quadratic detrend of the same record.
        q = csv_rows(run_cli("analyze", df, "--remove-drift", "quadratic", "--devs", "adev", "--no-ci").out)
        ref = adev(detrend(read_phase(df); method=:quadratic), [1, 2, 4, 8, 16, 32])
        @test isapprox(parse.(Float64, getindex.(q, 5))[1:6], ref.dev; rtol=1e-6)
        r = run_cli("check", df, "--remove-drift", "linear", "--format", "json")
        @test occursin("\"removed_drift_per_day\":null", r.out)
        @test occursin("\"removed_offset\":null", run_cli("check", df, "--format", "json").out)
        @test run_cli("check", df, "--remove-drift", "cubic").code == 1
        # Gaps (NaN) are skipped by the fit, not spread through it.
        xg = copy(X); xg[100] = NaN
        gf = write_table(joinpath(dir, "gap.csv"), xg .+ 3e-11 .* t)
        r = run_cli("check", gf, "--remove-drift", "linear", "--format", "json")
        @test r.code == 0
        @test isapprox(parse(Float64, match(r"\"removed_offset\":([^,}]+)", r.out)[1]), 3e-11; rtol=0.1)
        # With --fillgaps, the slots put back by the gap filler stay out of
        # the fit: it equals a least-squares fit over the measured rows only.
        tg = setdiff(0:length(X)-1, 120:179)
        xr = X[tg .+ 1] .+ 2e-12 .* tg .^ 2
        hf = joinpath(dir, "holes.csv")
        open(hf, "w") do io
            for (ti, v) in zip(tg, xr); println(io, ti, ",", v); end
        end
        r = run_cli("check", hf, "--fillgaps", "--remove-drift", "quadratic", "--format", "json")
        @test r.code == 0
        c = (tg .- (length(X) - 1) / 2)
        coef = [ones(length(c)) c c .^ 2] \ xr
        @test isapprox(parse(Float64, match(r"\"removed_offset\":([^,}]+)", r.out)[1]), coef[2]; rtol=1e-6)
        @test isapprox(parse(Float64, match(r"\"removed_drift_per_day\":([^,}]+)", r.out)[1]), 2 * coef[3] * 86400; rtol=1e-6)
        # One column with NaN readings, filled: the filled readings stay out of the fit too.
        nf = joinpath(dir, "nan.txt")
        xn = X .+ 3e-11 .* (0:length(X)-1); xn[100:139] .= NaN
        write(nf, join(xn, "\n"))
        r = run_cli("check", nf, "--tau0", "1", "--fillgaps", "--remove-drift", "linear", "--format", "json")
        ok = isfinite.(xn); c = (0:length(X)-1)[ok] .- (length(X) - 1) / 2
        @test isapprox(parse(Float64, match(r"\"removed_offset\":([^,}]+)", r.out)[1]),
                       ([ones(length(c)) c] \ xn[ok])[2]; rtol=1e-6)
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

    @testset "live" begin
        floats(s) = [parse(Float64, v) for v in split(s, ',') if !isempty(v)]
        readings(out) = reduce(vcat, [floats(m[1]) for m in eachmatch(r"\"readings\":\[([^\]]*)\]", out)])
        adevs(out) = floats(match(r"\"adev\",\"dev\":\[([^\]]*)\]", split(strip(out), '\n')[end])[1])

        r = run_cli("live", "sim", "--fast", "--count", "64", "--every", "32", "--m", "1,2", "--devs", "adev")
        @test r.code == 0
        lines = split(strip(r.out), '\n')
        @test length(lines) == 2
        @test occursin("\"n\":32,", lines[1]) && occursin("\"done\":false", lines[1])
        @test occursin("\"n\":64,", lines[2]) && occursin("\"done\":true", lines[2])
        @test adevs(r.out) ≈ adev(PhaseData(readings(r.out), 1.0), [1, 2]; ci=false).dev

        # A counter log replayed from a file: Hz readings against the first one.
        y = 1e-11 .* sinpi.((1:40) ./ 3.7)
        hf = joinpath(dir, "counter.txt")
        write(hf, join(string.(10e6 .* (1 .+ y)), "\n"))
        r = run_cli("live", hf, "--hz", "--tau0", "2", "--every", "100", "--m", "1,2,4")
        @test r.code == 0
        @test occursin("\"nominal\":1.0e7,", r.out)
        @test adevs(r.out) ≈ adev(FrequencyData(y, 2.0), [1, 2, 4]; ci=false).dev rtol = 1e-5

        # An LXI instrument: the SCPI query goes out before every reading. The
        # instrument runs in its own process, since the CLI blocks while it reads.
        instrument = open(`$(Base.julia_cmd()) --startup-file=no -e '
            using Sockets
            server = listen(Sockets.localhost, 0)
            println(Int(getsockname(server)[2])); flush(stdout)
            c = accept(server)
            while readline(c) == "MEAS:FREQ?"
                println(c, "+1.0000000001E+07")
            end'`)
        port = readline(instrument)
        r = run_cli("live", "tcp://127.0.0.1:" * port, "--query", "MEAS:FREQ?", "--hz",
                    "--nominal", "1e7", "--count", "5", "--m", "1")
        wait(instrument)
        @test r.code == 0
        @test readings(r.out) == fill(1.0000000001e7, 5)

        mf = joinpath(dir, "missed.txt")
        write(mf, "10000000.1\n9.91E+37\n10000000.2\n")
        r = run_cli("live", mf, "--hz")
        @test r.code == 2 && occursin("missed gate", r.err)

        @test run_cli("live").code == 1
        @test run_cli("live", "sim", "--freq").code == 1
        @test run_cli("live", "sim", "--devs", "totdev").code == 1
        @test run_cli("live", "sim", "--fillgaps").code == 1
        @test run_cli("live", "tcp://nowhere").code == 1
        r = run_cli("live", "tcp://127.0.0.1:" * port, "--count", "1")
        @test r.code == 2
        @test occursin("\"exit_code\":2}", r.out)
    end

    @testset "--range and τ limits" begin
        # Seconds 10 to 200 of a τ₀ = 1 s record are samples 11 to 201.
        rows = csv_rows(run_cli("analyze", pf, "--range", "10:200", "--devs", "adev", "--no-ci").out)
        ref = adev(PhaseData(X[11:201], 1.0), tau_values(Octave, 191, :adev); ci=false)
        @test parse.(Float64, getindex.(rows, 5)) == ref.dev
        r = run_cli("check", pf, "--range", "100:", "--format", "json")
        @test occursin("\"n\":$(length(X) - 100),", r.out) && occursin("\"range\":[100.0,null]", r.out)
        rows = csv_rows(run_cli("analyze", pf, "--devs", "adev", "--tau-min", "2", "--tau-max", "16").out)
        @test parse.(Float64, getindex.(rows, 4)) == [2.0, 4.0, 8.0, 16.0]
        @test run_cli("analyze", pf, "--range", "5").code == 1
        @test run_cli("analyze", pf, "--range", "9:3").code == 1
        @test run_cli("analyze", pf, "--range", "1e9:").code == 3
        @test run_cli("analyze", pf, "--tau-min", "16", "--tau-max", "2").code == 1
        @test run_cli("analyze", pf, "--devs", "adev", "--tau-min", "1.5", "--tau-max", "1.9").code == 3
        # The detrend is fitted to the range alone.
        rows = csv_rows(run_cli("analyze", pf, "--range", "10:200", "--detrend", "linear", "--devs", "adev", "--no-ci").out)
        ref = adev(detrend(PhaseData(X[11:201], 1.0); method=:linear), tau_values(Octave, 191, :adev); ci=false)
        @test parse.(Float64, getindex.(rows, 5)) ≈ ref.dev rtol=1e-12
        # Dynamic window centers stay in the whole record's time.
        whole = run_cli("dynamic", pf, "--window", "64", "--format", "json").out
        part = run_cli("dynamic", pf, "--window", "64", "--range", "64:", "--format", "json").out
        t(s) = parse(Float64, match(r"\"t\":\[([^,\]]+)", s)[1])
        @test t(part) == t(whole) + 64
        r = run_cli("report", pf, "--range", "10:200", "--tau-max", "16").out
        @test occursin("seconds 10 to 200", r) && occursin("up to 16 s", r)
    end

    @testset "report" begin
        r = run_cli("report", pf, "--devs", "adev,tdev", "--m", "1,2,4")
        @test r.code == 0 && startswith(r.out, "<!doctype html>") && endswith(r.out, "</html>\n")
        @test count("<svg", r.out) == 2                       # one plot per axis
        ref = adev(PhaseData(X, 1.0), [1, 2, 4])
        io = IOBuffer(); SigmaTauCLI.sci(io, ref.dev[2], 4)
        @test occursin(String(take!(io)), r.out)
        @test occursin("<pre>sigmatau report " * pf, r.out)
        @test run_cli("report", pf, "--format", "csv").code == 1
        sp = joinpath(dir, "it's here.csv"); cp(pf, sp; force=true)
        @test occursin("<pre>sigmatau report '" * replace(dirname(sp), "'" => "") * "/it'\\''s here.csv'", run_cli("report", sp).out)
    end

    @testset "synthetic record" begin
        r = run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev", "--format", "json")
        @test r.code == 0
        @test occursin("\"source\":\"synthetic\"", r.out)
        @test run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev").out ==
              run_cli("analyze", "--synth", "1024", "--seed", "7", "--devs", "adev").out
    end

    @testset "simulated clocks" begin
        # Clock i of --clocks K is column i of the ensemble the timescale uses.
        cc, _ = simulate_ensemble([ClockModel(; tau0=1.0, q1=1e-24 * 4.0^j, q2=1e-30 * 4.0^j) for j in 0:2], 512;
                                  reference=0, rng=Random.Xoshiro(9))
        for i in 1:3
            r = run_cli("check", "--synth", "512", "--clocks", "3", "--clock", string(i), "--seed", "9", "--series", "512", "--format", "json")
            @test r.code == 0
            m = match(r"\"series\":\{[^}]*?\"x\":\[([^\]]*)\]", r.out)
            @test m !== nothing
            @test parse.(Float64, split(m.captures[1], ','))[end] ≈ cc.x[end, i] rtol=1e-9
        end
        a = run_cli("analyze", "--synth", "1024", "--clocks", "3", "--clock", "1", "--devs", "adev").out
        b = run_cli("analyze", "--synth", "1024", "--clocks", "3", "--clock", "2", "--devs", "adev").out
        @test a != b                                        # the clocks differ
        @test run_cli("analyze", "--synth", "1024", "--clocks", "3", "--clock", "2", "--devs", "adev").out == b
        @test run_cli("check", "--synth", "64", "--clocks", "3", "--clock", "4").code == 1      # no such clock
        @test run_cli("check", "--synth", "64", "--clocks", "1").code == 1                      # needs two or more
        @test run_cli("check", "--synth", "64", "--clock", "2").code == 1                       # needs --clocks
        @test run_cli("check", pf, "--clocks", "3").code == 1                                   # synthetic only
    end

    @testset "noise models" begin
        ms = [ClockModel(; tau0=1.0, R=1e-22, q1=1e-24, q2=1e-30), ClockModel(; tau0=1.0, q1=4e-24)]
        cc, _ = simulate_ensemble(ms, 512; reference=0, rng=Random.Xoshiro(9))
        r = run_cli("check", "--synth", "512", "--seed", "9", "--noise-model", "1e-22:1e-24:1e-30:0", "--noise-model", "0:4e-24:0:0",
                    "--clock", "2", "--series", "512", "--format", "json")
        @test r.code == 0
        m = match(r"\"series\":\{[^}]*?\"x\":\[([^\]]*)\]", r.out)
        @test parse.(Float64, split(m.captures[1], ','))[end] ≈ cc.x[end, 2] rtol=1e-9
        @test run_cli("check", "--synth", "256", "--noise-model", "1e-22:1e-24:0:0").code == 0       # one clock
        # model: the closed form of the library
        r = run_cli("model", "--synth", "4096", "--noise-model", "1e-22:1e-24:1e-30:0", "--noise-model", "0:0:0:1e-40", "--format", "json")
        @test r.code == 0
        tau = [1.0, 2.0, 4.0]
        @test occursin("\"mdev\":[" * join(string.(powerlaw_sigma(ms[1], tau; deviation=:mdev)), ",")[1:12], r.out)
        @test occursin("\"mdev\":null,\"adev\":null,\"mhdev\":[", r.out)                       # random-run FM: Hadamard only
        @test run_cli("model", "--synth", "64").code == 1                                         # needs a model
        @test run_cli("model", pf).code == 1                                                      # and reads no file
        one = ["model", "--synth", "64", "--noise-model", "0:1e-24:0:0"]
        @test occursin("does not apply", run_cli(one..., "--fillgaps").err)
        mp = joinpath(dir, "model.json")
        @test isempty(run_cli(one..., "--out", mp, "--quiet").out)
        @test read(mp, String) == run_cli(one...).out
        @test run_cli("check", "--synth", "64", "--noise-model", "1:2:3").code == 1
        @test run_cli("check", "--synth", "64", "--noise-model", "-1:0:0:0").code == 1
        @test run_cli("check", "--synth", "64", "--clocks", "3", "--noise-model", "0:1e-24:0:0").code == 1
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
        # analyze still runs the gap-tolerant estimators and leaves the spectrum out
        r = run_cli("analyze", gf, "--devs", "adev", "--spectrum", "--format", "json")
        @test r.code == 0
        @test occursin("\"spectrum\":{\"error\":\"the spectrum needs a record without gaps\"}", r.out)
        @test occursin("\"deviation\":\"adev\"", r.out)
        # A one-column frequency file with NaN readings is analyzed once filled.
        yf = joinpath(dir, "y.txt")
        write(yf, join([i == 50 ? NaN : 1e-12 * sinpi(i / 7.3) for i in 1:256], "\n"))
        @test occursin("--fillgaps", run_cli("analyze", yf, "--freq", "--tau0", "1").err)
        @test run_cli("analyze", yf, "--freq", "--tau0", "1", "--devs", "adev", "--fillgaps").code == 0
        # A zero in a Stable32 frequency file is a gap, unless it is kept as a reading.
        zf = joinpath(dir, "z.frd")
        write(zf, join([i == 50 ? 0.0 : 1e-12 * sinpi(i / 7.3) for i in 1:256], "\n"))
        @test occursin("--fillgaps", run_cli("analyze", zf, "--tau0", "1").err)
        @test run_cli("analyze", zf, "--tau0", "1", "--devs", "adev", "--keep-zeros").code == 0
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
                         ["analyze", pf, "--devs", "totdev,mtotdev,theo1", "--bias", "stable32"],
                         ["analyze", pf, "--devs", "adev,htotdev,theoh", "--noise", "FLFM"],
                         ["analyze", pf, "--devs", "mdev,mhtotdev", "--noise", "fit"],
                         ["analyze", pf, "--format", "json", "--spectrum", "--series", "100"],
                         ["check", pf, "--format", "json"],
                         ["spectrum", pf],
                         ["cross", pf, "--with", pf, "--format", "json"],
                         ["analyze", joinpath(@__DIR__, "..", "..", "test", "fixtures", "formats", "timelab.tim")],
                         ["check", joinpath(@__DIR__, "..", "..", "test", "fixtures", "formats", "gpsdo.csv"), "--format", "json"],
                         ["dynamic", pf, "--dev", "dhdev"],
                         ["holdover", pf, "--train", "64", "--limit", "1e-9", "--fit", "--format", "json"],
                         ["holdover", pf, "--grow", "3", "--fit-drift", "--fit", "--format", "json"],
                         ["holdover", pf, "--train", "512", "--kalman", "--format", "json"],
                         ["holdover", joinpath(dir, "flicker.csv"), "--kalman", "--format", "json"],
                         ["steer", joinpath(dir, "steer.csv"), "--time-constant", "30", "--integral", "--format", "json", "--series", "40"],
                         ["steer", joinpath(dir, "steer.csv"), "--cost", "1:0:1e4", "--fit-drift"],
                         ["live", "sim", "--fast", "--count", "64", "--every", "32"],
                         ["analyze", pf, "--devs", "bogus", "--format", "json"],
                         ["model", "--synth", "64", "--noise-model", "1e-22:1e-24:1e-30:0"],
                         ["completion", "bash"],
                         ["help"],
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
