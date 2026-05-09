#!/usr/bin/env bash
# Fully-detached allantools synthetic-bench launcher.
set -euo pipefail
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/ian/SigmaTau.jl/.bench-venv/bin/python \
  /home/ian/SigmaTau.jl/benchmarks/bench/bench_allantools.py \
  --synth /home/ian/SigmaTau.jl/benchmarks/bench/synth \
  --out /home/ian/SigmaTau.jl/benchmarks/bench/results_allantools_synth.json \
  --m-max 512 \
  --kernels adev,mdev,hdev,tdev,totdev,mtotdev,htotdev
