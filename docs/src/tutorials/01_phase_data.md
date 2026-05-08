# Tutorial 1: Phase Data

How to construct a `PhaseData` record and load measurements from common
formats (DAT, CSV).

Narrative fills in a follow-up PR. Skeleton:

```@example phase
using SigmaTau
# build a synthetic phase record
x = randn(1000)
p = PhaseData(x, 1.0)
length(p.x), p.tau0
```
