# Paper and figures

`paper.md` and `paper.bib` are a draft software paper in the Journal of Open
Source Software format (not submitted). The `Paper` workflow renders it to PDF
with the JOSS toolchain whenever either file changes.

Figures about SigmaTau.jl, kept for papers and talks. They were drawn for
a poster, so they use a large type size, the Carlito
font and a crimson accent (`#9E1B32`). Data series use a colorblind-safe
palette (crimson, blue, amber, green, violet, validated for adjacent-pair
separation under protan, deutan and tritan simulation) with distinct marker
shapes and direct labels as a second channel.

| File | What it is |
|---|---|
| `figures/F1…F7` | the figures, PNG at 300 dpi and PDF |
| `figures/F8_explorer*.png` | screenshots of the GLMakie explorer, which has since been removed; kept as they are |
| `julia_figdata.jl` | runs the package itself to write `data/suite_fixture.csv` (every estimator on the validation fixture with noise type, EDF and interval), `data/dynamic_map.csv` (a DADEV map) and `data/streaming_check.csv` (streamed vs batch deviations) |
| `make_figures.py`, `make_schematic.py` | draw `figures/F1…F7` from those CSVs, `data/validation.json`, `bench/README.md` and `tools/artifacts/mhtotdev_mc_full.json` |
| `parse_validation.py` | parses the per-τ tables of `docs/src/explanation/validation.md` into `data/validation.json` |
| `stab_np.py`, `check_port.py` | a numpy port of ADEV/MDEV/HDEV/TOTDEV, checked against the Stable32 fixture values; a fallback for drawing Fig. 1 without Julia |

Every number in the figures comes from `README.md`, `docs/src/explanation/validation.md`,
`docs/src/mhtotdev.md`, `bench/README.md`, `tools/artifacts/mhtotdev_mc_full.json`,
`src/edf.jl` (`_coeff_mhtot`), `src/streaming.jl`, or the CSVs above.

Rebuild, from the repository root:

```bash
julia --project=. paper/julia_figdata.jl        # package output -> paper/data/*.csv
pip install numpy matplotlib
python3 paper/make_figures.py && python3 paper/make_schematic.py
```
