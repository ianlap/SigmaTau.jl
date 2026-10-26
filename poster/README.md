# Annual-review poster draft

A 40 in × 32 in, three-column poster about SigmaTau.jl for an expert
audience (industry and government partners in time and frequency). Alabama
crimson (`#9E1B32`) and white; data series use a colorblind-safe palette
(crimson, blue, amber, green, violet, validated for adjacent-pair separation
under protan, deutan and tritan simulation) with distinct marker shapes and
direct labels as a second channel.

| File | What it is |
|---|---|
| `SigmaTau_poster_draft.pptx` | the poster, one 40 × 32 in slide; text boxes, images and native tables, ready to move into the university template |
| `SigmaTau_poster_draft.pdf`, `SigmaTau_poster_preview.png` | renders for viewing |
| `content.js` | every heading, paragraph, bullet, caption, table cell and reference; edit here (`_{x}` subscript, `^{x}` superscript) |
| `build_poster.js` | lays `content.js` and `figures/` out on the slide and prints how full each column is |
| `content_v1_educational.js` | the first, background-heavy draft, kept for reference |
| `julia_figdata.jl` | runs the package itself to write `data/suite_fixture.csv` (every estimator on the validation fixture with noise type, EDF and interval), `data/dynamic_map.csv` (a DADEV map) and `data/streaming_check.csv` (streamed vs batch deviations) |
| `make_figures.py`, `make_schematic.py` | draw `figures/F1…F7` (PNG 300 dpi and PDF) from those CSVs, `data/validation.json`, `bench/README.md` and `tools/artifacts/mhtotdev_mc_full.json` |
| `figures/F8_explorer*.png` | screenshots written by the desktop explorer's headless smoke mode: `SIGMATAU_SMOKE_DIR=<dir> xvfb-run -a -s "-screen 0 1600x1050x24" julia --project=app app/StabilityExplorer.jl --smoke` |
| `parse_validation.py` | parses the per-τ tables of `docs/src/validation.md` into `data/validation.json` |
| `stab_np.py`, `check_port.py` | a numpy port of ADEV/MDEV/HDEV/TOTDEV, checked against the Stable32 fixture values; a fallback for drawing Fig. 1 without Julia |

Every number on the poster comes from `README.md`, `docs/src/validation.md`,
`docs/src/mhtotdev.md`, `bench/README.md`, `tools/artifacts/mhtotdev_mc_full.json`,
`src/edf.jl` (`_coeff_mhtot`), `src/streaming.jl`, or the CSVs above.

Rebuild, from the repository root:

```bash
julia --project=. poster/julia_figdata.jl        # package output -> poster/data/*.csv
pip install numpy matplotlib
python3 poster/make_figures.py && python3 poster/make_schematic.py
(cd poster && npm install && node build_poster.js)
```
