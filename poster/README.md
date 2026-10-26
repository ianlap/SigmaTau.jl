# Annual-review poster draft

A 40 in × 32 in, three-column poster about SigmaTau.jl, built from data that is
already in the repository. Alabama crimson (`#9E1B32`) and white; data series
use a colorblind-safe palette (crimson, blue, amber, green, violet, validated
for adjacent-pair separation under protan, deutan and tritan simulation) with
distinct marker shapes and direct labels as a second encoding.

| File | What it is |
|---|---|
| `SigmaTau_poster_draft.pptx` | the poster, one 40 × 32 in slide (drop the pieces into the university template) |
| `SigmaTau_poster_draft.pdf` | the same, rendered for viewing |
| `content.js` | every heading, paragraph, caption and table cell; edit here |
| `build_poster.js` | lays `content.js` and `figures/` out on the slide (`node build_poster.js`) |
| `figures/F1…F6` | PNG (300 dpi) and PDF of each figure |
| `make_figures.py`, `make_schematic.py` | regenerate the figures (`python3 make_figures.py`) |
| `stab_np.py`, `check_port.py` | numpy port of ADEV/MDEV/HDEV/TOTDEV used only to draw Fig. 1 when Julia is unavailable; `check_port.py` verifies it against the Stable32 fixture values (≤ 5 × 10⁻⁵, the fixture's printed precision) |
| `parse_validation.py`, `data/validation.json` | the per-τ tables of `docs/src/validation.md` as data |

Sources of every number on the poster: `README.md`, `docs/src/validation.md`,
`docs/src/mhtotdev.md`, `bench/README.md`, `tools/artifacts/mhtotdev_mc_full.json`,
`src/edf.jl` (`_coeff_mhtot`), `src/streaming.jl`.

Rebuild:

```bash
pip install numpy matplotlib            # figures
npm install pptxgenjs                   # poster (any directory on NODE_PATH, or edit the require path)
python3 poster/make_figures.py && python3 poster/make_schematic.py
node poster/build_poster.js
```
