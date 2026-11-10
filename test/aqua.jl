# test/aqua.jl — package hygiene checks that the General registry expects:
# method ambiguities, stale dependencies, compat bounds, piracy.
using Aqua, SigmaTau

Aqua.test_all(SigmaTau)
