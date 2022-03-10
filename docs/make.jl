using Documenter, Literate, Aeroelasticity

# Pre-install matplotlib
import Plots; Plots.pyplot()

# Generate examples
include("generate.jl")

# Build Documentation
makedocs(;
    modules = [Aeroelasticity],
    pages = [
        "Home" => "index.md",
        "Getting Started" => joinpath("examples", "guide.md"),
        "Examples" => [
            joinpath("examples", "section-stability.md"),
            joinpath("examples", "section-stability.md"),
        ],
        "Model Documentation" => [
            "Aerodynamic Models" => [
                "Steady Thin Airfoil Theory" => joinpath("models", "aerodynamics", "steady.md"),
                "Quasi-Steady Thin Airfoil Theory" => joinpath("models", "aerodynamics", "quasisteady.md"),
                "Wagner's Function" => joinpath("models", "aerodynamics", "wagner.md"),
                "Peters' Finite State" => joinpath("models", "aerodynamics", "peters.md"),
            ],
            "Structural Models" => [
                "Typical Section" => joinpath("models", "structures", "section.md"),
                "Geometrically Exact Beam Theory" => joinpath("models", "structures", "gxbeam.md"),
            ],
            "Coupled Models" => [
                "Steady Thin Airfoil Theory + Typical Section" => joinpath("models", "couplings", "steady-section.md"),
                "Quasi-Steady Thin Airfoil Theory + Typical Section" => joinpath("models", "couplings", "quasisteady-section.md"),
                "Wagner's Function + Typical Section" => joinpath("models", "couplings", "wagner-section.md"),
                "Peters' Finite State + Typical Section" => joinpath("models", "couplings", "peters-section.md"),
            ],
        ],
        "Developer Guide" => "developer.md",
        "API Reference" => [
            "Public" => joinpath("reference", "public.md"),
            "Internals" => joinpath("reference", "internals.md"),
        ],
    ],
    sitename = "Aeroelasticity.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

# Deploy Documentation
deploydocs(
    repo = "github.com/byuflowlab/Aeroelasticity.jl.git",
    devbranch = "main",
)
