using Documenter
using DocumenterCitations
using SigmaTau

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
            "Relativistic PNT" => [
                "theory/relativistic_clocks.md",
                "theory/relativistic_frames_and_timescales.md",
                "theory/relativistic_corrections.md",
                "theory/lunar_pnt_systems.md",
            ],
            "theory/ensembles_and_oscillator_networks.md",
        ],
        "Tutorials"       => [
            "tutorials/01_phase_data.md",
            "tutorials/02_compute_adev.md",
            "tutorials/03_identify_noise.md",
            "tutorials/04_confidence_intervals.md",
            "tutorials/05_single_clock_steering.md",
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
