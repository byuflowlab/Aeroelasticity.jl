module AerostructuralDynamics

using LinearAlgebra
using SpecialFunctions
using LinearMaps
using StaticArrays
import DifferentialEquations.ODEFunction

export AbstractModel, AerodynamicModel, StructuralModel, DynamicsModel
export TypicalSection
export PetersFiniteState

export number_of_states
export number_of_inputs
export number_of_parameters
export isinplace
export has_mass_matrix
export constant_mass_matrix
export linear_input_dependence
export defined_state_jacobian

export get_inputs
export get_mass_matrix
export get_state_jacobian

include("interface.jl")
include("structures/section.jl")
include("aerodynamics/peters.jl")
include("couplings/section.jl")

end
