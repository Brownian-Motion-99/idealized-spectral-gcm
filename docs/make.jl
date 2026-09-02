using Documenter, JGCM

DocMeta.setdocmeta!(JGCM, :DocTestSetup, :(using JGCM); recursive = true)

makedocs(;
    modules = [JGCM],
    doctest = true,
    linkcheck = false,
    checkdocs = :none,
    warnonly = [:missing_docs, :cross_references],
    authors = "Hao Chang <chunhaoc777@gmail.com>",
    repo = "https://github.com/Brownian-Motion-99/idealized-spectral-gcm.git",
    sitename = "idealized-spectral-gcm",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Brownian-Motion-99.github.io/idealized-spectral-gcm",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Model formulation" => [
            "Dynamical core" => "dynamics.md",
            "Physical parameterizations" => "physics.md",
        ],
        "User guide" => [
            "Model configuration" => "config.md",
            "Output and restarts" => "output.md",
        ],
        "Reference" => [
            "Notation" => "notation.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/Brownian-Motion-99/idealized-spectral-gcm.git",
    devbranch = "main",
    target = "build",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
)
