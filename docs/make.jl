using Documenter
using DocumenterCitations
using Literate
using SigmaTau
using Plots

Plots.gr()

# Render examples/*.jl into docs/src/tutorials/ (generated pages are gitignored).
const EXAMPLES_DIR  = joinpath(@__DIR__, "..", "examples")
const TUTORIALS_DIR = joinpath(@__DIR__, "src", "tutorials")
mkpath(TUTORIALS_DIR)
for jl in sort(readdir(EXAMPLES_DIR; join=true))
    endswith(jl, ".jl") || continue
    Literate.markdown(jl, TUTORIALS_DIR; documenter=true, credit=false)
end

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

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
        "Home" => "index.md",
        "Tutorials" => [
            "Getting started" => "getting_started.md",
            "tutorials/01_drifting_clock_hadamard.md",
            "tutorials/02_reading_your_data.md",
            "tutorials/03_three_cornered_hat.md",
        ],
        "How-to guides" => [
            "Choosing an estimator" => "guide/choosing.md",
            "Plotting and tables" => "guide/plotting.md",
            "Using the explorer" => "guide/explorer.md",
            "Command line" => "guide/cli.md",
        ],
        "Background" => [
            "Noise, EDF and confidence" => "explanation/noise_and_confidence.md",
            "Validation" => "explanation/validation.md",
            "MHTOTDEV bias and EDF" => "explanation/mhtotdev.md",
        ],
        "Reference" => [
            "Types" => "reference/types.md",
            "Estimators" => "reference/estimators.md",
            "Suites, τ grids and streaming" => "reference/suites.md",
            "Data and files" => "reference/data.md",
            "Preprocessing" => "reference/preprocessing.md",
            "Spectra" => "reference/spectra.md",
            "Noise and statistics" => "reference/statistics.md",
        ],
        "Bibliography" => "bibliography.md",
    ],
    doctest = true,
)

deploydocs(
    repo         = "github.com/ianlap/SigmaTau.jl.git",
    push_preview = true,
    devbranch    = "main",
)
