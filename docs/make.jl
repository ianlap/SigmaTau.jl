using Documenter
using SigmaTau, SigmaTauBase, SigmaTauStability, SigmaTauEnsemble

DocMeta.setdocmeta!(SigmaTauBase,      :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauStability, :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauEnsemble,  :DocTestSetup, :(using SigmaTau); recursive=true)

makedocs(
    sitename = "SigmaTau.jl",
    modules  = [SigmaTau, SigmaTauBase, SigmaTauStability, SigmaTauEnsemble],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://ianlap.github.io/SigmaTau.jl",
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Theory"          => [
            "theory/overview.md",
            "theory/allan_family.md",
            "theory/total_family.md",
            "theory/confidence.md",
            "theory/noise_id.md",
        ],
        "Tutorials"       => [
            "tutorials/01_phase_data.md",
            "tutorials/02_compute_adev.md",
            "tutorials/03_identify_noise.md",
            "tutorials/04_confidence_intervals.md",
            "tutorials/05_single_clock_steering.md",
        ],
        "API Reference"   => [
            "reference/base.md",
            "reference/stability.md",
            "reference/ensemble.md",
        ],
        "Validation"      => [
            "validation/methodology.md",
            "validation/stable32.md",
        ],
    ],
    doctest  = true,
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

deploydocs(
    repo         = "github.com/ianlap/SigmaTau.jl.git",
    push_preview = true,
    devbranch    = "main",
)
