"""
    AbstractModel

Supertype for all models.
"""
abstract type AbstractModel end

"""
    NoStateModel <: AbstractModel

Supertype for all models which contain no state variables.
"""
abstract type NoStateModel <: AbstractModel end

abstract type InPlaceness end
struct InPlace <: InPlaceness end
struct OutOfPlace <: InPlaceness end

# definition for models with no states
inplaceness(::Type{T}) where T<:NoStateModel = OutOfPlace()

# definition for combinations of models
function inplaceness(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isinplace(inplaceness(model_types...)) || any(isinplace.(model_types))
        return InPlace()
    else
        return OutOfPlace()
    end
end

isinplace(model::T) where T = isinplace(inplaceness(T))
isinplace(::Type{T}) where T = isinplace(inplaceness(T))
isinplace(::OutOfPlace) = false
isinplace(::InPlace) = true

abstract type MatrixType end
struct Empty <: MatrixType end
struct Zeros <: MatrixType end
struct Identity <: MatrixType end
struct Constant <: MatrixType end
struct Varying <: MatrixType end

# definitions for models with no states
mass_matrix_type(::Type{T}) where T<:NoStateModel = Empty()
state_jacobian_type(::Type{T}) where T<:NoStateModel = Empty()
input_jacobian_type(::Type{T}) where T<:NoStateModel = Empty()

# definitions for combinations of models
function mass_matrix_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isempty(mass_matrix_type(model_types...)) &&
        all(isempty.(mass_matrix_type.(model_types)))
        return Empty()
    elseif iszero(mass_matrix_type(model_types...)) &&
        all(iszero.(mass_matrix_type.(model_types)))
        return Zeros()
    elseif iszero(mass_matrix_type(model_types...)) &&
        all(isidentity.(mass_matrix_type.(model_types)))
        return Identity()
    elseif isconstant(mass_matrix_type(model_types...)) &&
        all(isconstant.(input_jacobian_type.(model_types))) &&
        all(isconstant.(mass_matrix_type.(model_types)))
        return Constant()
    else
        return Varying()
    end
end

function state_jacobian_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isempty(state_jacobian_type(model_types...)) &&
        all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif iszero(state_jacobian_type(model_types...)) &&
        all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif iszero(state_jacobian_type(model_types...)) &&
        all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif isconstant(state_jacobian_type(model_types...)) &&
        all(isconstant.(input_jacobian_type.(model_types))) &&
        all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    else
        return Varying()
    end
end

function input_jacobian_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if all(isempty.(input_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(input_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(input_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(input_jacobian_type.(model_types)))
        return Constant()
    else
        return Varying()
    end
end

isempty(::MatrixType) = false
isempty(::Empty) = true

iszero(::MatrixType) = false
iszero(::Zeros) = true

isidentity(::MatrixType) = false
isidentity(::Identity) = true

isconstant(::MatrixType) = false
isconstant(::Empty) = true
isconstant(::Zeros) = true
isconstant(::Constant) = true

abstract type InputDependence end
struct Linear <: InputDependence end
struct Nonlinear <: InputDependence end

input_dependence_type(::Type{T}) where T<:AbstractModel = Nonlinear()
input_dependence_type(::Type{T}) where T<:NoStateModel = Linear()

function linear_input_dependence(model)
    return _linear_input_dependence(input_dependence_type(typeof(model)))
end
_linear_input_dependence(::Linear) = true
_linear_input_dependence(::Nonlinear) = false
