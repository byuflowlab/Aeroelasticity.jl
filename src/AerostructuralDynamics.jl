module AerostructuralDynamics

using LinearAlgebra
using LinearMaps
using StaticArrays
import DifferentialEquations.ODEFunction

export AbstractModel
export TypicalSection
export Steady, QuasiSteady, Wagner, Peters

export number_of_states
export number_of_inputs

export number_of_parameters
export isinplace
export has_mass_matrix
export constant_mass_matrix
export linear_input_dependence
export defined_state_jacobian

export state_indices
export input_indices
export parameter_indices

export get_inputs, get_inputs!
export get_mass_matrix, get_mass_matrix!
export get_rates, get_rates!
export get_state_jacobian, get_state_jacobian!

include("traits.jl")
include("interface.jl")
include("structures/section.jl")
include("aerodynamics/quasisteady.jl")
include("aerodynamics/wagner.jl")
include("aerodynamics/peters.jl")

end
