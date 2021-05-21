using Documenter, AerostructuralDynamics

makedocs(;
    modules = [AerostructuralDynamics],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Developer Guide" => "developer.md",
        "Library" => "library.md",
    ],
    sitename = "AerostructuralDynamics.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
    # format = LaTeX(), # uncomment for PDF output
)

deploydocs(
    repo = "github.com/byuflowlab/AerostructuralDynamics.jl.git",
)
