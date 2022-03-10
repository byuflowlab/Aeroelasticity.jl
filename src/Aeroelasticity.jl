module Aeroelasticity

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
export number_of_states, number_of_inputs, number_of_parameters
export state_indices, input_indices, parameter_indices
export separate_states, separate_parameters
export linearize, get_eigen, correlate_eigenmodes

export Steady, QuasiSteady, Wagner, Peters, LiftingLine
export Section, GXBeamAssembly

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
include(joinpath("models", "aerodynamics", "liftingline.jl"))

include(joinpath("models", "structures", "section.jl"))
include(joinpath("models", "structures", "gxbeam.jl"))

include(joinpath("models", "couplings", "steady-section.jl"))
include(joinpath("models", "couplings", "quasisteady-section.jl"))
include(joinpath("models", "couplings", "wagner-section.jl"))
include(joinpath("models", "couplings", "peters-section.jl"))

end
