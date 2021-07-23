module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using StaticArrays
using DiffEqBase

# Models
export AbstractModel
export NoStateModel

export Steady
export QuasiSteady
export Wagner
export Peters
export LiftingLine
export TypicalSection
export RigidBody
export GEBT
export couple_models

# Interface
export number_of_states
export number_of_inputs
export number_of_parameters
export state_indices
export input_indices
export parameter_indices
export get_inputs, get_inputs!
export get_rates, get_rates!
export get_mass_matrix, get_mass_matrix!
export get_state_jacobian, get_state_jacobian!
export get_eigen
export get_ode

# Helper Functions
export set_states
export set_inputs
export set_parameters

# Developer Interface
export AbstractModel
export NoStateModel
export inplaceness
export mass_matrix_type
export state_jacobian_type
export input_jacobian_type
export get_input_jacobian
export get_input_mass_matrix
export get_input_mass_matrix!
export get_input_state_jacobian
export get_input_state_jacobian!

# needed for traits
import Base.isempty
import Base.iszero

# Interface Functions
include("traits.jl")
include("interface.jl")

# 2D Aerodynamic Models
include("aerodynamics/quasisteady.jl")
include("aerodynamics/wagner.jl")
include("aerodynamics/peters.jl")

# 3D Aerodynamic Models
include("aerodynamics/liftingline.jl")

# 2D Structural Models
include("structures/section.jl")

# 3D Structural Models
include("structures/rigidbody.jl")
include("structures/gxbeam.jl")

# 2D Coupled Models
include("couplings/quasisteady-section.jl")
include("couplings/quasisteady-liftingline.jl")
include("couplings/wagner-section.jl")
include("couplings/wagner-liftingline.jl")
include("couplings/peters-section.jl")
include("couplings/peters-liftingline.jl")

# 3D Coupled Models
include("couplings/liftingline-rigidbody.jl")
include("couplings/liftingline-gxbeam.jl")
include("couplings/liftingline-gxbeam-rigidbody.jl")

end
