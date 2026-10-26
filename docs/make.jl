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
        "Home"       => "index.md",
        "Plotting and tables" => "plotting.md",
        "Validation" => "validation.md",
        "MHTOTDEV bias and EDF" => "mhtotdev.md",
        "Tutorials"  => [
            "tutorials/01_drifting_clock_hadamard.md",
            "tutorials/02_reading_your_data.md",
            "tutorials/03_three_cornered_hat.md",
        ],
        "API" => [
            "reference/types.md",
            "reference/stab.md",
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
