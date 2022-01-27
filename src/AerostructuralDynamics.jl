module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using FillArrays
using StaticArrays
using SciMLBase
using Printf

import Base.isempty
import Base.iszero

export assemble_model
export assemble_states, assemble_states!
export assemble_parameters, assemble_parameters!
export separate_states, separate_parameters
export linearize, get_eigen, correlate_eigenmodes

export Steady, QuasiSteady, Wagner, Peters, LiftingLine
export Section

export Submodel, Coupling

export section_coordinates

include("util.jl")
include("models.jl")
include("jacobians.jl")
include("internals.jl")
include("interface.jl")

include(joinpath("models", "aerodynamics", "steady.jl"))
include(joinpath("models", "aerodynamics", "quasisteady.jl"))
include(joinpath("models", "aerodynamics", "wagner.jl"))
include(joinpath("models", "aerodynamics", "peters.jl"))
include(joinpath("models", "structures", "section.jl"))
include(joinpath("models", "couplings", "steady-section.jl"))
include(joinpath("models", "couplings", "quasisteady-section.jl"))
include(joinpath("models", "couplings", "wagner-section.jl"))
include(joinpath("models", "couplings", "peters-section.jl"))

end
