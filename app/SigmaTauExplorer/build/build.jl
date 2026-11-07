# Build the desktop explorer as a standalone, relocatable app folder with
# PackageCompiler.create_app. No Julia is needed on the machine that runs it.
#
#   julia app/SigmaTauExplorer/build/build.jl [outdir]
#
# Default outdir: app/SigmaTauExplorer/dist/SigmaTauExplorer. The launcher is
# bin/sigmatau-explorer (bin\sigmatau-explorer.exe on Windows). Takes tens of
# minutes and needs a display (a desktop session, or xvfb-run on Linux).
# Build on each target OS; PackageCompiler does not cross-compile.
using Pkg

const PKG = normpath(joinpath(@__DIR__, ".."))
const OUT = abspath(get(ARGS, 1, joinpath(PKG, "dist", "SigmaTauExplorer")))

# Resolve the app's own environment (SigmaTau comes from this checkout).
Pkg.activate(PKG)
Pkg.instantiate()

# PackageCompiler lives in a throwaway environment, not in the app.
Pkg.activate(; temp=true)
Pkg.add(name="PackageCompiler", version="2")
using PackageCompiler

t = @elapsed create_app(PKG, OUT;
    executables = ["sigmatau-explorer" => "julia_main"],
    precompile_execution_file = joinpath(@__DIR__, "precompile.jl"),
    incremental = false,
    # Ship only the stdlibs the dependency graph declares. If the app dies
    # with a missing stdlib, set this to false and file the missing dep.
    filter_stdlibs = true,
    include_lazy_artifacts = true,
    force = true)

bytes = sum(filesize(joinpath(r, f)) for (r, _, fs) in walkdir(OUT) for f in fs)
println("Built ", OUT, " in ", round(t / 60; digits=1), " min, ",
        round(bytes / 2^20; digits=1), " MiB on disk")
