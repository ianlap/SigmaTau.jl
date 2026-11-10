# Validation fixtures

| Path | What | Source |
|---|---|---|
| `stable32gen.DAT`, `stable32out/`, `allantools_out/` | 8192-point record with Stable32 and allantools runs | used by `test/stab/` |
| `noise/alpha_*.txt` | power-law phase records, α = +2 … −4, N = 4096, τ0 = 1 s | `tools/validation/make_noise_fixtures.jl` (Kasdin–Walter, fixed seeds) |
| `stable32_phasedat/PHASE.DAT` | Stable32's PHASE.DAT (the NBS 1000-point frequency set integrated to phase) | distributed with Stable32; copied from the allantools test suite |
| `stable32_phasedat/phase_dat_*.txt` | Stable32 1.53 / 1.60 result files on PHASE.DAT (σ, α, 68.3 % bounds) | allantools test suite (`tests/phasedat/`) |
| `stable32_flickerpm/flickerpm_8192` | 8192-point flicker PM frequency record, τ0 = 1 s, as Stable32 wrote it | Stable32 noise generator, flicker PM coefficient 1, 2026-11-09 |
| `stable32_flickerpm/flickerpm_8192_tables.txt` | Stable32's printed tables on it: overlapping Allan, modified Allan, overlapping Hadamard, total, modified total, Hadamard total (octave AF 1–2048, total to 4096; 0.683 confidence) | Stable32; σ asserted by `test/validation/`, α and bounds not |
| `allantools/*.csv` | allantools 2024.6 and scipy runs on all of the above | `tools/validation/gen_allantools.py` |

The records are committed so both programs read identical numbers. The
Cross-check workflow regenerates `allantools/` with a live install and runs
`test/validation/` against the fresh tables.
