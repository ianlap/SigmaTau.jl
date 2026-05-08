# Validation Methodology

SigmaTau.jl is cross-validated against three independent references:

1. **Stable32** (W. Riley) — desktop application, the de facto industry
   reference for time-and-frequency stability analysis.
2. **allantools** (A. Wallin) — Python library, second independent
   numerical reference.
3. **AllanLab** (MATLAB) — third reference, locked-in fixture.

Three-way agreement defines the rtol floor; documented disagreements are
boundary-policy differences (TOTDEV/HTOTDEV/MTOTDEV reflection conventions).

Detailed comparison narrative lands in a follow-up PR.
