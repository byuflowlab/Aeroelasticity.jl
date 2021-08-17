module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using StaticArrays
using DiffEqBase
using RecipesBase
using Printf

# Models
export Steady
export QuasiSteady
export Wagner
export Peters
export LiftingLine
export TypicalSection
export RigidBody
export GEBT
export LinearFlap
export LiftingLineFlaps
export Trim

# Interface
export couple_models
export number_of_states
export number_of_inputs
export number_of_parameters
export state_indices
export input_indices
export parameter_indices
export get_states, set_states!
export get_inputs, set_inputs!
export get_parameters, set_parameters!
export get_additional_parameters, set_additional_parameters!
export separate_states
export separate_inputs
export separate_parameters
export get_inputs, get_inputs!
export get_rates, get_rates!
export get_mass_matrix, get_mass_matrix!
export get_state_jacobian, get_state_jacobian!
export get_eigen
export get_ode
export correlate_eigenmodes

# needed for traits
import Base.isempty
import Base.iszero

# Interface Functions
include("traits.jl")
include("interface.jl")
include("mode-tracking.jl")

# Two-Dimensional Models

# Aerodynamic Models
include(joinpath("aerodynamics", "quasisteady.jl"))
include(joinpath("aerodynamics", "wagner.jl"))
include(joinpath("aerodynamics", "peters.jl"))

# Structural Models
include(joinpath("structures", "section.jl"))
include(joinpath("structures", "liftingline-section.jl"))

# Control Surface Models
include(joinpath("control-surfaces", "linearflap.jl"))

# Controller Models
include(joinpath("controllers", "liftingline-section-control.jl"))

# Coupled Models
include(joinpath("couplings", "two-dimensional", "quasisteady-section.jl"))
include(joinpath("couplings", "two-dimensional", "wagner-section.jl"))
include(joinpath("couplings", "two-dimensional", "peters-section.jl"))

include(joinpath("couplings", "two-dimensional", "quasisteady-section-linearflap.jl"))
include(joinpath("couplings", "two-dimensional", "wagner-section-linearflap.jl"))
include(joinpath("couplings", "two-dimensional", "peters-section-linearflap.jl"))

include(joinpath("couplings", "two-dimensional", "quasisteady-liftingline.jl"))
include(joinpath("couplings", "two-dimensional", "wagner-liftingline.jl"))
include(joinpath("couplings", "two-dimensional", "peters-liftingline.jl"))

include(joinpath("couplings", "two-dimensional", "quasisteady-liftingline-linearflap.jl"))
include(joinpath("couplings", "two-dimensional", "wagner-liftingline-linearflap.jl"))
include(joinpath("couplings", "two-dimensional", "peters-liftingline-linearflap.jl"))

# Three-Dimensional Models

# Aerodynamic Models
include(joinpath("aerodynamics", "liftingline.jl"))

# Structural Models
include(joinpath("structures", "gxbeam.jl"))

# Dynamics Models
include(joinpath("dynamics", "rigidbody.jl"))

# Control Surface Models
include(joinpath("control-surfaces", "liftinglineflaps.jl"))

# Controller Models
include(joinpath("controllers", "trim.jl"))

# Coupled Models
include(joinpath("couplings", "three-dimensional", "liftingline-rigidbody.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam-rigidbody.jl"))

include(joinpath("couplings", "three-dimensional", "liftingline-rigidbody-liftinglineflaps.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam-liftinglineflaps.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam-rigidbody-liftinglineflaps.jl"))

include(joinpath("couplings", "three-dimensional", "liftingline-rigidbody-liftinglineflaps-trim.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam-liftinglineflaps-trim.jl"))
include(joinpath("couplings", "three-dimensional", "liftingline-gxbeam-rigidbody-liftinglineflaps-trim.jl"))

end
