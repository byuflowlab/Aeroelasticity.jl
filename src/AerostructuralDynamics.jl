module AerostructuralDynamics

using ArnoldiMethod
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using StaticArrays
using DiffEqBase

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

# Interface
export couple_models
export number_of_states
export number_of_inputs
export number_of_parameters
export state_indices
export input_indices
export parameter_indices
export set_states
export set_inputs
export set_parameters
export get_inputs, get_inputs!
export get_rates, get_rates!
export get_mass_matrix, get_mass_matrix!
export get_state_jacobian, get_state_jacobian!
export get_eigen
export get_ode

# needed for traits
import Base.isempty
import Base.iszero

# Interface Functions
include("traits.jl")
include("interface.jl")

# 2D Aerodynamic Models
include(joinpath("aerodynamics", "quasisteady.jl"))
include(joinpath("aerodynamics", "wagner.jl"))
include(joinpath("aerodynamics", "peters.jl"))

# 3D Aerodynamic Models
include(joinpath("aerodynamics", "liftingline.jl"))

# 2D Structural Models
include(joinpath("structures", "section.jl"))

# 3D Structural Models
include(joinpath("structures", "rigidbody.jl"))
include(joinpath("structures", "gxbeam.jl"))

# 2D Control Surface Models
include(joinpath("control-surfaces", "linearflap.jl"))

# 3D Control Surface Models
include(joinpath("control-surfaces", "liftinglineflaps.jl"))

# 2D Coupled Models
include(joinpath("couplings", "quasisteady-section.jl"))
include(joinpath("couplings", "quasisteady-section-linearflap.jl"))
include(joinpath("couplings", "wagner-section.jl"))
include(joinpath("couplings", "wagner-section-linearflap.jl"))
include(joinpath("couplings", "peters-section.jl"))
include(joinpath("couplings", "peters-section-linearflap.jl"))

# 2D to 3D Coupled Models (internal)
include(joinpath("couplings", "quasisteady-liftingline.jl"))
include(joinpath("couplings", "quasisteady-liftingline-linearflap.jl"))
include(joinpath("couplings", "wagner-liftingline.jl"))
include(joinpath("couplings", "wagner-liftingline-linearflap.jl"))
include(joinpath("couplings", "peters-liftingline.jl"))
include(joinpath("couplings", "peters-liftingline-linearflap.jl"))

# 3D Coupled Models
include(joinpath("couplings", "liftingline-rigidbody.jl"))
include(joinpath("couplings", "liftingline-rigidbody-liftinglineflaps.jl"))
include(joinpath("couplings", "liftingline-gxbeam.jl"))
include(joinpath("couplings", "liftingline-gxbeam-liftinglineflaps.jl"))
include(joinpath("couplings", "liftingline-gxbeam-rigidbody.jl"))
include(joinpath("couplings", "liftingline-gxbeam-rigidbody-liftinglineflaps.jl"))

end
