# SigmaTauExplorer: the desktop explorer as a standalone app

This package wraps `app/StabilityExplorer.jl` so it can be compiled into an
app that runs without a Julia install. The explorer code is not duplicated
here; `src/SigmaTauExplorer.jl` includes it.

## Build an app folder (PackageCompiler)

```bash
julia app/SigmaTauExplorer/build/build.jl                 # Linux: prefix with xvfb-run -a on a headless box
app/SigmaTauExplorer/dist/SigmaTauExplorer/bin/sigmatau-explorer
```

The build traces `build/precompile.jl` (the `--smoke` self-test plus a hidden
window) into a custom sysimage, so the window opens without compiling. It
needs a real display while building: on a headless CI runner GLMakie's
OpenGL paths are only partly traced, and the app starts slower. Build on each
target OS; there is no cross-compilation. The launcher takes the same options
as the script (`--fontsize=N`, `--size=WxH`, `--smoke`). For `--smoke` outside
a source checkout, set `SIGMATAU_SMOKE_FIXTURE` to a copy of
`test/fixtures/validation/stable32gen.DAT`.

The manual GitHub Actions workflow `Explorer app` runs this build on Linux
under xvfb, runs the built app's `--smoke`, and uploads the folder as a
tarball with its size and timings in the job summary.

## Installers (AppBundler.jl)

The module also defines `@main`, so it runs as
`julia --project=app/SigmaTauExplorer -m SigmaTauExplorer`, which is the
launch convention [AppBundler.jl](https://github.com/PeaceFounder/AppBundler.jl)
packages into a DMG (macOS), MSIX (Windows) or Snap (Linux):

```bash
julia -e 'using Pkg; Pkg.Apps.add("AppBundler")'   # once; puts `appbundler` in ~/.julia/bin
appbundler build app/SigmaTauExplorer --build-dir=app/SigmaTauExplorer/dist --selfsign
```

Its sysimage mode and signing options are configured per its own docs
(<https://peacefounder.org/AppBundler.jl/>); that step has not been tried on
this package yet.
