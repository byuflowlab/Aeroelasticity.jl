module Aeroelasticity

using ArnoldiMethod
using Arpack
using FiniteDiff
using FLOWMath
using ForwardDiff
using GXBeam
using LinearAlgebra
using LinearMaps
using FillArrays
using IterativeSolvers
using StaticArrays
using SciMLBase
using SparseArrays
using SparseDiffTools
using Symbolics
using Printf
using UnPack

export CoupledModel

export number_of_states, state_indices, separate_states

export residual, residual!, state_jacobian, state_jacobian!, rate_jacobian, rate_jacobian!

export linearize, dense_eigen, sparse_eigen, correlate_eigenmodes

export Steady
export QuasiSteady
export Wagner
export Peters
export LiftingLine, LiftingLineParameters
export Section, section_coordinates
export RigidBody
export GXBeamAssembly, GXBeamInputs, GXBeamParameters
export LiftingLineGXBeamAssembly, LiftingLineGXBeamParameters
export LiftingLineRigidBody
export SteadySection, QuasiSteadySection, WagnerSection, PetersSection

include("models.jl")

include("interface.jl")

# aerodynamic models
include(joinpath("models", "aerodynamics", "steady.jl"))
include(joinpath("models", "aerodynamics", "quasisteady.jl"))
include(joinpath("models", "aerodynamics", "wagner.jl"))
include(joinpath("models", "aerodynamics", "peters.jl"))
include(joinpath("models", "aerodynamics", "liftingline.jl"))

# structural models
include(joinpath("models", "structures", "section.jl"))
include(joinpath("models", "structures", "liftingline.jl"))
include(joinpath("models", "structures", "rigidbody.jl"))
include(joinpath("models", "structures", "gxbeam.jl"))

# typical section couplings
include(joinpath("models", "couplings", "steady-section.jl"))
include(joinpath("models", "couplings", "quasisteady-section.jl"))
include(joinpath("models", "couplings", "wagner-section.jl"))
include(joinpath("models", "couplings", "peters-section.jl"))

# lifting line section couplings
include(joinpath("models", "couplings", "steady-liftingline.jl"))
include(joinpath("models", "couplings", "quasisteady-liftingline.jl"))
include(joinpath("models", "couplings", "wagner-liftingline.jl"))
include(joinpath("models", "couplings", "peters-liftingline.jl"))

# rigid body couplings
include(joinpath("models", "couplings", "liftingline-rigidbody.jl"))

# geometrically exact beam theory couplings
include(joinpath("models", "couplings", "liftingline-gxbeam.jl"))

end
