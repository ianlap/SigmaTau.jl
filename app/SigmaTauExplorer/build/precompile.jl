# Workload that PackageCompiler traces while building the sysimage: the same
# headless self-test as `StabilityExplorer.jl --smoke` (every page, every
# export, the file round trips), plus opening and closing a hidden window so
# the interactive screen path is compiled too. Needs a display: run the build
# on a desktop session or under xvfb-run.
using SigmaTauExplorer

withenv("SIGMATAU_SMOKE_DIR" => mktempdir()) do
    SigmaTauExplorer.explorer_main(["--smoke"])
end

try
    app = SigmaTauExplorer.build()
    app.reload(false)
    app.run()
    screen = display(app.fig; visible=false)
    sleep(2)
    close(screen)
catch err
    @warn "hidden-window warm-up failed; the first window will compile at launch" exception = err
end
