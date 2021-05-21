module AerostructuralDynamics

using LinearAlgebra
using SpecialFunctions
using LinearMaps
using StaticArrays
import DifferentialEquations.ODEFunction

export AbstractModel, AerodynamicModel, StructuralModel, DynamicsModel
export TypicalSection
export PetersFiniteState

include("interface.jl")
include("structures/section.jl")
include("aerodynamics/peters.jl")

end
