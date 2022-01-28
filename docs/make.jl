using Documenter, Aeroelasticity

makedocs(;
    modules = [Aeroelasticity],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Model Documentation" => [
            "Aerodynamic Models" => [
                "Steady Thin Airfoil Theory" => joinpath("aerodynamics", "steady.md"),
                "Quasi-Steady Thin Airfoil Theory" => joinpath("aerodynamics", "quasisteady.md"),
                "Wagner's Function" => joinpath("aerodynamics", "wagner.md"),
                "Peters' Finite State" => joinpath("aerodynamics", "peters.md"),
                "Lifting Line" => joinpath("aerodynamics", "liftingline.md"),
            ],
            "Structural Models" => [
                "Typical Section" => joinpath("structures", "section.md"),
                "Rigid Body" => joinpath("structures", "rigidbody.md"),
                "Geometrically Exact Beam Theory" => joinpath("structures", "gxbeam.md"),
            ],
            "Coupled Models" => [
                "Steady Thin Airfoil Theory + Typical Section" => joinpath("couplings", "steady-section.md"),
                "Quasi-Steady Thin Airfoil Theory + Typical Section" => joinpath("couplings", "quasisteady-section.md"),
                "Wagner's Function + Typical Section" => joinpath("couplings", "wagner-section.md"),
                "Peters' Finite State + Typical Section" => joinpath("couplings", "peters-section.md"),
                "Lifting Line + Rigid Body" => joinpath("couplings", "liftingline-rigidbody.md"),
                "Lifting Line + Geometrically Exact Beam Theory" => joinpath("couplings", "liftingline-gxbeam.md"),
                "Lifting Line + Geometrically Exact Beam Theory + Rigid Body" => joinpath("couplings", "liftingline-gxbeam-rigidbody.md"),
            ],
        ],
        "Library" => [
            "Public" => joinpath("library", "public.md"),
            "Internals" => joinpath("library", "internals.md"),
        ],
        "Developer Guide" => "developer.md",
    ],
    sitename = "Aeroelasticity.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
    # format = LaTeX(), # uncomment for PDF output
)

deploydocs(
    repo = "github.com/byuflowlab/Aeroelasticity.jl.git",
    devbranch = "main",
)
