module Aeroelasticity

using ArnoldiMethod
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

export assemble_model
export assemble_states, assemble_states!
export assemble_parameters, assemble_parameters!
export number_of_states, number_of_inputs, number_of_parameters
export state_indices, input_indices, parameter_indices
export separate_states, separate_parameters

export CoupledModel

export residual, residual!, state_jacobian, state_jacobian!, rate_jacobian, rate_jacobian!

export linearize, dense_eigen, sparse_eigen, correlate_eigenmodes

export Steady, QuasiSteady, Wagner, Peters, LiftingLine
export Section, RigidBody, GXBeamAssembly

export section_coordinates

include("models.jl")
include("interface.jl")

include(joinpath("models", "aerodynamics", "steady.jl"))
include(joinpath("models", "aerodynamics", "quasisteady.jl"))
include(joinpath("models", "aerodynamics", "wagner.jl"))
include(joinpath("models", "aerodynamics", "peters.jl"))
include(joinpath("models", "aerodynamics", "liftingline.jl"))

include(joinpath("models", "structures", "section.jl"))
include(joinpath("models", "structures", "rigidbody.jl"))
include(joinpath("models", "structures", "liftingline.jl"))
include(joinpath("models", "structures", "gxbeam.jl"))

include(joinpath("models", "couplings", "steady-section.jl"))
include(joinpath("models", "couplings", "quasisteady-section.jl"))
include(joinpath("models", "couplings", "wagner-section.jl"))
include(joinpath("models", "couplings", "peters-section.jl"))

include(joinpath("models", "couplings", "steady-liftingline.jl"))
include(joinpath("models", "couplings", "quasisteady-liftingline.jl"))
include(joinpath("models", "couplings", "wagner-liftingline.jl"))
include(joinpath("models", "couplings", "peters-liftingline.jl"))

include(joinpath("models", "couplings", "liftingline-gxbeam.jl"))

end
