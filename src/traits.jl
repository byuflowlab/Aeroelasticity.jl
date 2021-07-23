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

# --- Trait Types --- #

abstract type InPlaceness end
struct InPlace <: InPlaceness end
struct OutOfPlace <: InPlaceness end

abstract type MatrixType end
struct Empty <: MatrixType end
struct Zeros <: MatrixType end
struct Identity <: MatrixType end
struct Constant <: MatrixType end
struct Linear <: MatrixType end
struct Nonlinear <: MatrixType end

# --- Trait Functions --- #

"""
    inplaceness(::Type{T})

Return `InPlace()` if functions associated with model `T` are in-place
or `OutOfPlace()` if functions associated with model `T` are out-of-place.
"""
inplaceness(::Type{T}) where T

# models with no state variables use out-of-place definitions
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

"""
    inplaceness(::Type{T1}, ::Type{T2}, ..., ::Type{TN})

Return `InPlace()` if the functions associated with the input function
for coupled models `T1`, `T2`, ... `TN` are in-place or `OutOfPlace()`
if the functions associated with the input function for coupled models `T1`,
`T2`, ... `TN` are out-of-place.
"""
inplaceness(::Vararg{Type,N}) where N

"""
   mass_matrix_type(::Type{T})

Return
 - `Empty()`, if the mass matrix associated with model `T` is empty
 - `Zeros()`, if the mass matrix associated with model `T` is filled
    with zeros
 - `Identity()`, if the mass matrix associated with model `T` is the
    identity matrix
 - `Constant()`, if the mass matrix associated with model `T` is
    constant with respect to time
 - `Linear()`, if the mass matrix associated with model `T` may vary
    with respect to time

If no method is defined for the specified type, return `Linear()`.
"""
mass_matrix_type(::Type{T}) where T = Linear()

# models with no state variables have no mass matrix
mass_matrix_type(::Type{T}) where T<:NoStateModel = Empty()

# definition for combinations of models
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
        return Linear()
    end
end

"""
   state_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` is empty
 - `Zeros()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` is filled
    with zeros
 - `Identity()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` is the
    identity matrix
 - `Constant()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` is
    constant with respect to time
 - `Linear()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` may vary
    with respect to time, but is linear with respect to the states
 - `Nonlinear()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the state variables associated with model `T` may vary
    with respect to time, and is nonlinear with respect to the states

If no method is defined for the specified type, return `Nonlinear()`.
"""
state_jacobian_type(::Type{T}) where T =  Nonlinear()

# models with no state variables have no state jacobian
state_jacobian_type(::Type{T}) where T<:NoStateModel = Empty()

# definition for combinations of models
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
    elseif all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

"""
   input_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs is empty for model `T`
 - `Zeros()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs is filled with zeros for model `T`
 - `Identity()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs is the identity matrix for model `T`
 - `Constant()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs is constant with respect to time for model `T`
 - `Linear()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs may vary with respect to time for model `T`, and
    is linear with respect to the inputs
 - `Nonlinear()`, if the jacobian of the mass matrix multiplied state rates
    with respect to the inputs may vary with respect to time for model `T`, and
    is nonlinear with respect to the inputs

If no method is defined for the specified type, return `Nonlinear()`.
"""
input_jacobian_type(::Type{T}) where T = Nonlinear()

# models with no state variables have no input jacobian
input_jacobian_type(::Type{T}) where T<:NoStateModel = Empty()

# definition for combinations of models
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
    elseif all(islinear.(input_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- dispatch functions --- #

isinplace(model::T) where T = isinplace(inplaceness(T))
isinplace(::Type{T}) where T = isinplace(inplaceness(T))
isinplace(::OutOfPlace) = false
isinplace(::InPlace) = true

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

islinear(::MatrixType) = false
islinear(::Empty) = true
islinear(::Zeros) = true
islinear(::Constant) = true
islinear(::Linear) = true
