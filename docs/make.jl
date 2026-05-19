using Documenter
using DocumenterCitations
using Literate
using SigmaTau
# Plot backend: PGFPlotsX renders LaTeX-quality vector PDFs and
# font-matches the docs body. Loaded here so that any `@example` block
# that subsequently `using Plots` picks up PGFPlotsX as the default
# backend automatically (Plots.jl module state is process-global).
# Requires `pdflatex` / `lualatex` and `pdftocairo` in PATH — see
# .github/workflows/Documentation.yml for the CI install of texlive
# packages and poppler-utils.
using Plots
using PGFPlotsX
Plots.pgfplotsx()
# Enable `\text{…}` (and friends) inside math labels — PGFPlotsX's
# default preamble ships pgfplots only, not amsmath.
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amsmath}")

# ── Literate.jl: render examples/*.jl into docs/src/tutorials/ ──────────────
# Each top-level `examples/*.jl` is single-source: edit the script,
# rebuild, and the matching `tutorials/<name>.md` regenerates. The
# generated pages are gitignored (only the `.jl` files are tracked).
const EXAMPLES_DIR  = joinpath(@__DIR__, "..", "examples")
const TUTORIALS_DIR = joinpath(@__DIR__, "src", "tutorials")
mkpath(TUTORIALS_DIR)
for jl in sort(readdir(EXAMPLES_DIR; join=true))
    endswith(jl, ".jl") || continue
    Literate.markdown(jl, TUTORIALS_DIR;
                      documenter = true,
                      credit     = false)
end

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :authoryear,
)

DocMeta.setdocmeta!(SigmaTau, :DocTestSetup, :(using SigmaTau); recursive=true)

makedocs(
    sitename = "SigmaTau.jl",
    modules  = [SigmaTau, SigmaTau.Stab, SigmaTau.Est],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://ianlap.github.io/SigmaTau.jl",
        mathengine = Documenter.MathJax3(),
    ),
    plugins = [bib],
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Theory"          => [
            "theory/overview.md",
            "theory/allan_family.md",
            "theory/total_family.md",
            "theory/confidence.md",
            "theory/noise_id.md",
            "theory/ensemble_overview.md",
            "theory/kalman.md",
            "theory/steering.md",
            "theory/clock_ensembles.md",
            "theory/validation.md",
        ],
        "Tutorials"       => [
            "tutorials/00_julia_for_metrologists.md",
            "tutorials/01_phase_data.md",
            "tutorials/02_compute_adev.md",
            "tutorials/03_kalman_single_clock.md",
            "tutorials/04_kalman_pid_steering.md",
            "tutorials/05_holdover_comparison.md",
            "tutorials/06_three_cornered_hat.md",
            "tutorials/07_clock_ensemble.md",
        ],
        "API Reference"   => [
            "reference/types.md",
            "reference/stab.md",
            "reference/est.md",
        ],
        "Validation"      => [
            "validation/methodology.md",
            "validation/stable32.md",
        ],
        "Bibliography"    => "bibliography.md",
    ],
    doctest  = true,
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

deploydocs(
    repo         = "github.com/ianlap/SigmaTau.jl.git",
    push_preview = true,
    devbranch    = "main",
)
