---
title: 'SigmaTau.jl: validated time-and-frequency stability analysis in Julia'
tags:
  - Julia
  - frequency stability
  - Allan deviation
  - atomic clocks
  - time metrology
authors:
  - name: Ian Lapinski
    affiliation: 1
affiliations:
  - name: The University of Alabama, United States
    index: 1
date: 8 November 2026
bibliography: paper.bib
---

# Summary

Clocks, oscillators and frequency references are characterized by how their
fractional frequency wanders as a function of averaging time τ. The standard
statistics for this are the Allan variance [@allan-1966-statistics] and its
modified, Hadamard, total and Thêo relatives, defined in IEEE Std 1139
[@ieee1139-2022-definitions] and collected in NIST SP 1065
[@riley-2008-sp1065]. A stability estimate is only useful with its
uncertainty, which depends on the dominant power-law noise at each τ and on
the equivalent degrees of freedom (EDF) of the estimator
[@greenhall-2003-edf-stability].

SigmaTau.jl is a Julia [@bezanson-2017-julia] package that computes these
statistics end to end. One call on a phase or frequency record returns, at
every τ, the deviation, the identified noise type
[@riley-2004-lag1-acf], the EDF from the Greenhall–Riley algorithm and the χ²
confidence interval. It covers the Allan (`adev`, `mdev`, `tdev`), Hadamard
(`hdev`, `mhdev`, `htdev`), total (`totdev`, `mtotdev`, `ttotdev`, `htotdev`,
`mhtotdev`) [@greenhall-1999-totvar; @howe-1999-modtotvar;
@howe-2005-tothvar-ieee] and Thêo families, time-error statistics (`mtie`,
`tierms`, `pdev`), dynamic deviations, streaming accumulators
[@dobrogowski-2007-realtime-adev], N-cornered-hat and Groslambert-covariance
separation with confidence intervals [@vernotte-2019-tch-gcov-uncertainty;
@lantz-2019-klts], cross-spectra [@baudiquez-2020-cross-spectrum-statistics], spectral
estimates and power-law noise synthesis. The same core is compiled into a
stand-alone command line, `sigmatau`, which starts in about 30 ms, and drives a
desktop explorer and a Python wrapper.

# Statement of need

Two tools dominate the field. Stable32 [@riley-stable32] is the reference
application, but it is a closed-source Windows GUI that cannot be scripted
into an analysis pipeline. allantools [@wallin-allantools] is an open Python
library that reports deviations at full precision, but leaves noise
identification, EDF and confidence intervals to the user, applies no bias
correction to the total estimators, and runs the O(N·m) total-family kernels
slowly enough that five days of 1 s data takes half an hour. In Julia,
AllanDeviations.jl provides the basic deviations without uncertainty.

SigmaTau targets the laboratory and operations user who needs Stable32's
answers, including its confidence intervals and bias conventions, from code:
in a notebook, in CI on a clock ensemble, or in a compiled binary. Every
default is chosen to match Stable32, and every place where the defaults
differ from a reference is a named keyword that reproduces it.

# Validation

Agreement with the references is the package's central claim, so it is
tested, not asserted in prose. The validation test group compares SigmaTau
row by row against Stable32 output, allantools output, the NBS 1000-point
reference values [@riley-2008-sp1065] and scipy's Welch spectra; 7037
comparisons are asserted in CI, and a separate CI job regenerates the
allantools references from source. On the shipped 8192-point fixture the core
estimators agree with Stable32 to its five printed significant figures
(≤ 5·10⁻⁵ relative) and with allantools to ≤ 10⁻¹³; the χ² bounds agree with
Stable32 to ≤ 0.05 % wherever the two programs identify the same noise type.
The comparison is published as a generated report, with each residual
traced to a definitional choice (a bias-correction policy, an EDF table edge
or a noise-identification convention) rather than tolerated.

Where no reference exists the package documents its own: MHTOTDEV has no
published bias factor or EDF, so both are measured by a reproducible Monte
Carlo harness against synthesized power-law noise and shipped with the
methodology.

# Performance

The kernels are written for Julia's compiler and multithreaded. On identical
inputs `adev`, `mdev`, `hdev` and `tdev` run 10 to 40 times faster than
allantools, `totdev` about 70 times, and `mtotdev`/`htotdev` 1,600 to 4,200
times faster (0.46 s against 30 min on a 407,000-sample record). The
benchmark harness and both sets of measurements are in the repository.

# Acknowledgements

<!-- Funding and acknowledgements to be added by the author. -->

# References
