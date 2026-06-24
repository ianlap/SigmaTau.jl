#!/usr/bin/env julia
# bin/sigmatau.jl — non-interactive command-line launcher for SigmaTau.jl.
#
# Run from the package directory:
#   julia --project=. bin/sigmatau.jl --help
#   julia --project=. bin/sigmatau.jl analyze data.DAT --phase --tau0 1 \
#         --devs adev,mdev,hdev --tau octave --out results.csv
#
# All logic lives in `SigmaTau.cli_main`, which returns a process exit code; this
# wrapper just forwards ARGS and exits with that code. The function is designed
# to be PackageCompiler-friendly so this same entry point can become a compiled
# binary later (out of scope for now).
import SigmaTau
exit(SigmaTau.cli_main(ARGS))
