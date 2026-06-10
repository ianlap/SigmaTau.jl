# Validation Methodology

SigmaTau.jl is cross-validated against independent references where those
references exist:

1. **Stable32** (W. Riley) — desktop application, the de facto industry
   reference for time-and-frequency stability analysis.
2. **allantools** (A. Wallin) — Python library, second independent
   numerical reference.

Agreement with the shipped fixtures defines the rtol floor; documented
disagreements are boundary-policy differences (TOTDEV/HTOTDEV/MTOTDEV
reflection conventions). Estimators without an external implementation are
checked against inlined legacy kernels and internal identities.

Validation is about *numerical agreement*, not speed. For head-to-head
timings against allantools — including the ~4,000× speedup on the
modified-total kernels for long records — see [Performance](../performance.md).

Detailed comparison narrative lands in a follow-up PR.
