using Documenter, JGCM

DocMeta.setdocmeta!(JGCM, :DocTestSetup, :(using JGCM); recursive=true)

const page_rename = Dict("developer.md" => "Developer docs")

function nice_name(file)
    file = replace(file, r"^[0-9]*-" => "")
    if haskey(page_rename, file)
        return page_rename[file]
    end
    return splitext(file)[1] |> x -> replace(x, "-" => " ") |> titlecase
end

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
    pages=[
        "Home" => "index.md",
        [
            nice_name(file) => file for
            file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
        ]... # Added splatting (...) to ensure the array is flattened into the pages list
    ],
)

deploydocs(; 
    repo = "github.com/Brownian-Motion-99/idealized-spectral-gcm.git", 
    devbranch = "v0.2.0",
    target = "build",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)