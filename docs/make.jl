using Documenter
using DocumenterCitations
using Literate
using SigmaTau
using Plots
using PGFPlotsX

# Plot backend: PGFPlotsX renders LaTeX-quality vector PDFs and font-matches the
# docs body when a TeX engine is available. Local builds without TeX fall back to
# GR so docs can still be checked without installing the full LaTeX toolchain.
# CI installs lualatex and pdftocairo in .github/workflows/Documentation.yml.
if any(cmd -> Sys.which(cmd) !== nothing, ("lualatex", "pdflatex", "xelatex"))
    Plots.pgfplotsx()
    # Enable `\text{…}` (and friends) inside math labels — PGFPlotsX's default
    # preamble ships pgfplots only, not amsmath.
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amsmath}")
else
    @warn "No LaTeX engine found; using GR for local docs plot rendering"
    Plots.gr()
end

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
    modules  = [SigmaTau],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://ianlap.github.io/SigmaTau.jl",
        mathengine = Documenter.MathJax3(),
    ),
    plugins = [bib],
    checkdocs = :exports,
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "From Stable32"   => "migration_from_stable32.md",
        "From allantools" => "migration_from_allantools.md",
        "Performance"     => "performance.md",
        "Plotting"        => "plotting.md",
        "FAQ"             => "faq.md",
        "Theory"          => [
            "theory/overview.md",
            "theory/allan_family.md",
            "theory/total_family.md",
            "theory/confidence.md",
            "theory/noise_id.md",
            "theory/spectral.md",
            "theory/mhtotdev_bias_edf.md",
            "theory/pdev_confidence.md",
            "theory/validation.md",
        ],
        "Tutorials"       => [
            "tutorials/00_julia_for_metrologists.md",
            "tutorials/01_phase_data.md",
            "tutorials/02_compute_adev.md",
            "tutorials/03_drifting_clock_hadamard.md",
            "tutorials/04_reading_your_data.md",
            "tutorials/06_three_cornered_hat.md",
        ],
        "API Reference"   => [
            "reference/types.md",
            "reference/stab.md",
        ],
        "Validation"      => [
            "validation/methodology.md",
            "validation/stable32.md",
        ],
        "Bibliography"    => "bibliography.md",
    ],
    doctest  = true,
)

deploydocs(
    repo         = "github.com/ianlap/SigmaTau.jl.git",
    push_preview = true,
    devbranch    = "main",
)
