# Tools

Scripts that produce numbers the package and its documentation rest on.
The tests do not run them.

## Exact bias and EDF models

`mhtot_models.jl` computes the exact bias and equivalent degrees of freedom
of MHTOTDEV on Gaussian power-law noise and fits the EDF model. It is the
source of `_MHTOT_B`, `_MHTOT_B_SMALL`, `_coeff_mhtot` and
`_MHTOT_WHPM_EDF_SMALL` in `src/edf.jl`. No published table covers this
estimator, and
[MHTOTDEV bias and EDF](https://ianlap.github.io/SigmaTau.jl/dev/explanation/mhtotdev/)
describes the method.

```bash
julia -t auto tools/mhtot_models.jl       # about three minutes on four cores
```

`total_pm_models.jl` does the same for TOTDEV and HTOTDEV on white and
flicker PM, which the published tables do not cover, and for the bias of
MTOTDEV and HTOTDEV below m = 16. It is the source of the PM rows of
`_total_bias`, of `_totvar_pm_edf` and `_htot_pm_edf`, and of
`_MTOT_B_SMALL` in `src/edf.jl`
([Total estimators on PM noise](https://ianlap.github.io/SigmaTau.jl/dev/explanation/total_pm/)).

```bash
julia -t auto tools/total_pm_models.jl    # under a minute
```

Both scripts build the quadratic form of each statistic from its
definition, apart from the kernels of the package, and use no random
numbers. `test/stab/mhtotdev_models.jl` and `test/stab/total_pm_models.jl`
hold the kernels to the numbers they print.

## `validation/`

The scripts behind the cross-check of `test/validation/`
([Validation](https://ianlap.github.io/SigmaTau.jl/dev/explanation/validation/)).
The tables they write are committed under `test/fixtures/validation/`, so
they are run only to regenerate them:

```bash
# the power-law records
julia --project=. tools/validation/make_noise_fixtures.jl

# the tables of allantools and scipy, about a minute
pip install allantools==2024.6 scipy numpy
python3 tools/validation/gen_allantools.py

# the table of AllanDeviations.jl
julia --project=tools/validation/allandeviations -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/validation/allandeviations tools/validation/gen_allandeviations.jl

# REPORT.md, which lists every comparison, and the tables of the Validation page
julia --project=. tools/validation/report.jl
julia --project=. tools/validation/tables.jl
```

`gen_allantools.py` and `gen_allandeviations.jl` write the committed tables
again byte for byte (run with allantools 2024.06, numpy 2.5.4, scipy 1.18.1
and AllanDeviations.jl 0.3.0). With `--out DIR`, and with a directory as
the argument of the Julia script, they write elsewhere, which is how the
Cross-check workflow compares fresh tables with the committed ones.

`make_noise_fixtures.jl` need not reproduce its records: on another machine
a few samples of each record came out different in the last printed digit,
since the synthesis does not round alike everywhere. The tables of the
other libraries are computed from those records, so write them again after
it.
