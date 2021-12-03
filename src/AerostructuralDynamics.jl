module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using FillArrays
using StaticArrays
using DiffEqBase
using RecipesBase
using Printf

# Interface
export number_of_states, number_of_inputs, number_of_parameters
export state_indices, input_indices, parameter_indices
export get_states, get_inputs, get_parameters, get_additional_parameters
export set_states!, set_inputs!, set_parameters!, set_additional_parameters!
export separate_states, separate_inputs, separate_parameters
export get_residual, get_residual!
export get_rate_jacobian, get_rate_jacobian!
export get_state_jacobian, get_state_jacobian!
export get_input_jacobian, get_input_jacobian!
export get_parameter_jacobian, get_parameter_jacobian!
export get_time_gradient, get_time_gradient!
export get_coupling_inputs, get_coupling_inputs!
export linearize, get_eigen, get_ode, get_dae
export correlate_eigenmodes

# needed for traits
import Base.isempty
import Base.iszero

# Aerodynamic Models
export steady_model, quasisteady_model, wagner_model, peters_model, liftingline_model

# Structural Models
export typical_section_model, liftingline_section_model, gxbeam_model

# Dynamics Models
export rigidbody_model

# Coupled Models
export steady_section_model, quasisteady_section_model, wagner_section_model, peters_section_model
export steady_liftingline_model, quasisteady_liftingline_model, wagner_liftingline_model, peters_liftingline_model
export liftingline_rigidbody_model, liftingline_gxbeam_model, liftingline_gxbeam_rigidbody_model

# Helper Functions
include("util.jl")

# Model Types
include("models.jl")

# Traits
include("traits.jl")

# Interface Functions
include("interface.jl")

# Model Definitions
include(joinpath("models", "aerodynamics", "steady.jl"))
include(joinpath("models", "aerodynamics", "quasisteady.jl"))
include(joinpath("models", "aerodynamics", "wagner.jl"))
include(joinpath("models", "aerodynamics", "peters.jl"))
include(joinpath("models", "aerodynamics", "liftingline.jl"))
include(joinpath("models", "structures", "section.jl"))
include(joinpath("models", "structures", "liftingline.jl"))
include(joinpath("models", "structures", "gxbeam.jl"))
include(joinpath("models", "dynamics", "rigidbody.jl"))
include(joinpath("models", "coupled", "two-dimensional", "steady-section.jl"))
include(joinpath("models", "coupled", "two-dimensional", "quasisteady-section.jl"))
include(joinpath("models", "coupled", "two-dimensional", "wagner-section.jl"))
include(joinpath("models", "coupled", "two-dimensional", "peters-section.jl"))
include(joinpath("models", "coupled", "two-dimensional", "steady-liftingline.jl"))
include(joinpath("models", "coupled", "two-dimensional", "quasisteady-liftingline.jl"))
include(joinpath("models", "coupled", "two-dimensional", "wagner-liftingline.jl"))
include(joinpath("models", "coupled", "two-dimensional", "peters-liftingline.jl"))
include(joinpath("models", "coupled", "three-dimensional", "liftingline-gxbeam.jl"))
include(joinpath("models", "coupled", "three-dimensional", "liftingline-rigidbody.jl"))

end
