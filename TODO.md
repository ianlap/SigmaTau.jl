# SigmaTau.jl — Roadmap

Working list of outstanding engineering work, sorted by priority. Items move
from this file to `CHANGELOG.md` once landed.

> **Audit date**: 2026-05-07. Source-truth audit performed against the actual
> code (not stale logs). Several items previously listed as "blocking" in the
> superseded `next_step_prompts.md` were already implemented by a prior agent
> session; those are noted under "✅ Recently completed (verify)".

---

## 🟡 High (completeness)

---

## 🟡 Medium (correctness gaps)

- [ ] **FFT-based FLPM/FLFM noise synthesis** for the multi-noise validation
  testset. Currently `MTOTDEV across noise regimes` exercises WPM, WHFM, and
  RWFM (synthesizable without an FFT); FLPM and FLFM need a `1/f` filter
  (see `legacy/julia/src/noise_gen.jl`).
- [ ] **Published `_coeff_totvar` / `_coeff_htot` for α=2,1** if the
  ADEV/HDEV-style fallback turns out to disagree with Stable32 at WPM/FLPM.
  The fallback in [`stats/edf.jl`](lib/SigmaTauStability/src/stats/edf.jl)
  is the pragmatic substitute used in absence of a published value; replace
  with a real entry from a peer-reviewed source if/when one surfaces.

---

## 🟢 Low (polish)

- [ ] **`RelativisticClock`** — empty struct in
  [`clocks.jl`](lib/SigmaTauEnsemble/src/models/clocks.jl#L20). Lunar PNT
  future work; no concrete need yet.
- [ ] **`UDFactorizedFilter`, `KuramotoOscillator`** — empty structs in
  [`filters.jl`](lib/SigmaTauEnsemble/src/estimators/filters.jl). Stubs
  reserved for low-observability lunar-distance and SWaP-constrained
  nearest-neighbor estimators respectively.
- [ ] **PID steering port** — legacy `kalman_filter` includes a PID
  controller; the new `predict!`/`update!` loop is steering-free. Required
  before clock-steering examples can be ported.
- [ ] **More `examples/`** — only `examples/quickstart.jl` exists; add a
  Kalman-only example, a `FrequencyData`-vs-`PhaseData` walkthrough, and a
  multi-clock ensemble scenario when the multi-clock model lands.

---

## 🟢 Documentation

- [ ] Author `README.md` (project intro, install, quickstart, link to docs).
- [ ] Maintain `CHANGELOG.md` (Keep-a-Changelog format).
- [ ] Author user-facing docs (Documenter.jl) in `docs/`.

---

## ✅ Recently completed (since 2026-05-07 session start)

See [CHANGELOG.md](CHANGELOG.md) once authored. In-flight items are tracked in
this file; once a checkbox here gets ticked and the change is committed, move
the entry to the changelog and delete it from here.
