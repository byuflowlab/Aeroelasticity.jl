module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using StaticArrays
using DiffEqBase

export AbstractModel
export TypicalSection, LiftingLineSection, RigidBody, GEBT
export Steady, QuasiSteady, Wagner, Peters, LiftingLine

export couple_models
export stability_analysis
export ode_function

export number_of_states
export number_of_inputs

export number_of_parameters
export isinplace
export has_mass_matrix
export constant_mass_matrix
export defined_state_jacobian

export state_indices
export input_indices
export parameter_indices

export get_inputs, get_inputs!
export get_mass_matrix, get_mass_matrix!
export get_rates, get_rates!
export get_state_jacobian, get_state_jacobian!
export get_input_jacobian, get_input_jacobian!

export default_parameters, default_inputs

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
