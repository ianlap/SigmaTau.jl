# tools/

`mhtot_models.jl` computes the exact bias and EDF of MHTOTDEV and fits the EDF
model: the source of `_MHTOT_B` and `_coeff_mhtot` in `src/edf.jl`. MHTOTDEV
has no published table, so this is the provenance of those numbers; the method
is described in `docs/src/explanation/mhtotdev.md`.

```bash
julia -t auto tools/mhtot_models.jl    # about five minutes
```

## validation/

The reference cross-check suite behind `test/validation/`:

```bash
julia --project=. tools/validation/make_noise_fixtures.jl   # power-law records (once; committed)
pip install allantools==2024.6 scipy numpy
python3 tools/validation/gen_allantools.py                   # allantools / scipy tables (~5 min)
julia --project=. tools/validation/report.jl                 # REPORT.md: agreement and discrepancies
julia --project=. tools/validation/tables.jl                 # per-τ tables of the docs Validation page
```
