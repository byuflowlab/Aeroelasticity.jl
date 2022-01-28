using Aeroelasticity
using LinearAlgebra
using GXBeam
using ForwardDiff
using Test

include("testutil.jl")

include(joinpath("models", "aerodynamics", "steady.jl"))
include(joinpath("models", "aerodynamics", "quasisteady.jl"))
include(joinpath("models", "aerodynamics", "wagner.jl"))
include(joinpath("models", "aerodynamics", "peters.jl"))

include(joinpath("models", "structures", "section.jl"))

include(joinpath("models", "couplings", "steady-section.jl"))
include(joinpath("models", "couplings", "quasisteady-section.jl"))
include(joinpath("models", "couplings", "wagner-section.jl"))
include(joinpath("models", "couplings", "peters-section.jl"))
