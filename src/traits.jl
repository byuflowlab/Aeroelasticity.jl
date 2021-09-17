# --- Models --- #

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

"""
    InterfaceModel <: AbstractModel

Supertype for all models which are used to extend 2D models to 3D models.
"""
abstract type InterfaceModel <: AbstractModel end

"""
    ResidualModel <: AbstractModel

Supertype for all models whose state variables may be defined by solving a set
of residual equations.
"""
abstract type ResidualModel <: AbstractModel end

"""
    UnsteadyModel <: AbstractModel

Supertype for all models which represent dynamic systems of equations.
"""
abstract type UnsteadyModel <: AbstractModel end

# --- Inplaceness Trait --- #

abstract type InPlaceness end
struct InPlace <: InPlaceness end
struct OutOfPlace <: InPlaceness end

"""
    inplaceness(::Type{T})

Return `InPlace()` if functions associated with model `T` are in-place
or `OutOfPlace()` if functions associated with model `T` are out-of-place.
"""
inplaceness(::Type{T}) where T

# default definition for models with no state variables
inplaceness(::Type{T}) where T<:NoStateModel = OutOfPlace()

# default definition for coupled models
function inplaceness(::Type{T}) where T <: NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isinplace(coupling_inplaceness(model_types...)) || any(isinplace.(model_types))
        return InPlace()
    else
        return OutOfPlace()
    end
end

"""
    coupling_inplaceness(::Type{T1}, ::Type{T2}, ..., ::Type{TN})

Return `InPlace()` if the functions associated with the coupling function
for coupled models `T1`, `T2`, ... `TN` are in-place or `OutOfPlace()`
if the functions associated with the coupling function for coupled models `T1`,
`T2`, ... `TN` are out-of-place.
"""
coupling_inplaceness(::Vararg{Type,N}) where N

isinplace(model::T) where T = isinplace(inplaceness(T))
isinplace(::Type{T}) where T = isinplace(inplaceness(T))
isinplace(::InPlace) = true
isinplace(::OutOfPlace) = false

# --- Matrix Type Trait --- #

abstract type MatrixType end
struct Empty <: MatrixType end
struct Zeros <: MatrixType end
struct Identity <: MatrixType end
struct Constant <: MatrixType end
struct Invariant <: MatrixType end
struct Linear <: MatrixType end
struct Nonlinear <: MatrixType end

"""
    rate_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state rates is an empty matrix
 - `Zeros()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state rates is a zero matrix
 - `Identity()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state rates is the identity matrix
 - `Invariant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state rates is independent of the state rates,
    state variables, inputs, parameters, and time.
 - `Constant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state rates is independent of the state rates,
    state variables, inputs, and time.
 - `Linear()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state rates may vary with respect to time, and is
    linear with respect to the state rates.
 - `Nonlinear()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state rates may vary with respect to time,
    and is nonlinear with respect to the state rates.

If no method is defined for the specified type, return `Nonlinear()`.
"""
rate_jacobian_type(::Type{T}) where T = Nonlinear()

# models with no state variables have no rate jacobian
rate_jacobian_type(::Type{T}) where T <: NoStateModel = Empty()

# definition for combinations of models
function rate_jacobian_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isempty(coupling_rate_jacobian_type(model_types...)) &&
        all(isempty.(rate_jacobian_type.(model_types)))
        return Empty()
    elseif iszero(coupling_rate_jacobian_type(model_types...)) &&
        all(iszero.(rate_jacobian_type.(model_types)))
        return Zeros()
    elseif iszero(coupling_rate_jacobian_type(model_types...)) &&
        all(isidentity.(rate_jacobian_type.(model_types)))
        return Identity()
    elseif isinvariant(coupling_rate_jacobian_type(model_types...)) &&
        all(isinvariant.(input_jacobian_type.(model_types))) &&
        all(isinvariant.(rate_jacobian_type.(model_types)))
        return Invariant()
    elseif isconstant(coupling_rate_jacobian_type(model_types...)) &&
        all(isconstant.(input_jacobian_type.(model_types))) &&
        all(isconstant.(rate_jacobian_type.(model_types)))
        return Constant()
    else
        return Linear()
    end
end

"""
    state_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state variables is an empty matrix
 - `Zeros()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state variables is a zero matrix
 - `Identity()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state variables is the identity matrix
 - `Invariant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state variables is independent of the state rates,
    state variables, inputs, parameters, and time.
 - `Constant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state variables is independent of the state rates,
    state variables, inputs, and time.
 - `Linear()`, if the jacobian of the residual expression associated with model `T`
    with respect to the state variables may vary with respect to time, and is
    linear with respect to the states.
 - `Nonlinear()`, if the jacobian of the residual expression associated with
    model `T` with respect to the state variables may vary with respect to time,
    and is nonlinear with respect to the states.

If no method is defined for the specified type, return `Nonlinear()`.
"""
state_jacobian_type(::Type{T}) where T = Nonlinear()

# models with no state variables have no state jacobian
state_jacobian_type(::Type{T}) where T <: NoStateModel = Empty()

# definition for combinations of models
function state_jacobian_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if isempty(coupling_state_jacobian_type(model_types...)) &&
        all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif iszero(coupling_state_jacobian_type(model_types...)) &&
        all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif iszero(coupling_state_jacobian_type(model_types...)) &&
        all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif isinvariant(coupling_state_jacobian_type(model_types...)) &&
        all(isinvariant.(input_jacobian_type.(model_types))) &&
        all(isinvariant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif isconstant(coupling_state_jacobian_type(model_types...)) &&
        all(isconstant.(input_jacobian_type.(model_types))) &&
        all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif islinear(coupling_state_jacobian_type(model_types...)) &&
        all(islinear.(input_jacobian_type.(model_types))) &&
        all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

"""
    input_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the residual expression associated with model `T`
    with respect to the inputs is an empty matrix
 - `Zeros()`, if the jacobian of the residual expression associated with model `T`
    with respect to the inputs is a zero matrix
 - `Identity()`, if the jacobian of the residual expression associated with
    model `T` with respect to the inputs is the identity matrix
 - `Invariant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the inputs is independent of the state rates,
    state variables, inputs, parameters, and time.
 - `Constant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the inputs is independent of the state rates,
    state variables, inputs, and time.
 - `Linear()`, if the jacobian of the residual expression associated with model `T`
    with respect to the inputs may vary with respect to time, and is
    linear with respect to the inputs.
 - `Nonlinear()`, if the jacobian of the residual expression associated with
    model `T` with respect to the inputs may vary with respect to time,
    and is nonlinear with respect to the inputs.

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
    elseif all(isinvariant.(input_jacobian_type.(model_types)))
        return Invariant()
    elseif all(isconstant.(input_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(input_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

"""
    parameter_jacobian_type(::Type{T})

Return
 - `Empty()`, if the jacobian of the residual expression associated with model `T`
    with respect to the parameters is an empty matrix
 - `Zeros()`, if the jacobian of the residual expression associated with model `T`
    with respect to the parameters is a zero matrix
 - `Identity()`, if the jacobian of the residual expression associated with
    model `T` with respect to the parameters is the identity matrix
 - `Invariant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the parameters is independent of the state rates,
    state variables, inputs, parameters, and time.
 - `Constant()`, if the jacobian of the residual expression associated with
    model `T` with respect to the parameters is independent of the state rates,
    state variables, inputs, and time.
 - `Linear()`, if the jacobian of the residual expression associated with model `T`
    with respect to the parameters may vary with respect to time, and is
    linear with respect to the parameters.
 - `Nonlinear()`, if the jacobian of the residual expression associated with
    model `T` with respect to the parameters may vary with respect to time,
    and is nonlinear with respect to the parameters.

If no method is defined for the specified type, return `Nonlinear()`.
"""
parameter_jacobian_type(::Type{T}) where T = Nonlinear()

# models with no state variables have no input jacobian
parameter_jacobian_type(::Type{T}) where T <: NoStateModel = Empty()

# definition for combinations of models
function parameter_jacobian_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if all(isempty.(parameter_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(parameter_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(parameter_jacobian_type.(model_types)))
        return Identity()
    elseif all(isinvariant.(parameter_jacobian_type.(model_types)))
        return Invariant()
    elseif all(isconstant.(parameter_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(parameter_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

"""
    time_gradient_type(::Type{T})

Return
 - `Empty()`, if the derivative of the residual expression associated with model `T`
    with respect to time is an empty matrix
 - `Zeros()`, if the derivative of the residual expression associated with model `T`
    with respect to time is a zero matrix
 - `Invariant()`, if the derivative of the residual expression associated with
    model `T` with respect to time is independent of the state rates,
    state variables, inputs, parameters, and time.
 - `Constant()`, if the derivative of the residual expression associated with
    model `T` with respect to time is independent of the state rates,
    state variables, inputs, and time.
 - `Linear()`, if the derivative of the residual expression associated with model `T`
    with respect to time may vary with respect to time, and is
    linear with respect to time.
 - `Nonlinear()`, if the derivative of the residual expression associated with
    model `T` with respect to time may vary with respect to time,
    and is nonlinear with respect to time.

If no method is defined for the specified type, return `Nonlinear()`.
"""
time_gradient_type(::Type{T}) where T = Nonlinear()

# models with no state variables have no input jacobian
time_gradient_type(::Type{T}) where T <: NoStateModel = Empty()

# definition for combinations of models
function time_gradient_type(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    if all(isempty.(time_gradient_type.(model_types)))
        return Empty()
    elseif all(iszero.(time_gradient_type.(model_types)))
        return Zeros()
    elseif all(isinvariant.(time_gradient_type.(model_types)))
        return Invariant()
    elseif all(isconstant.(time_gradient_type.(model_types)))
        return Constant()
    elseif all(islinear.(time_gradient_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

coupling_rate_jacobian_type(model_types...) = Nonlinear()
coupling_state_jacobian_type(model_types...) = Nonlinear()
coupling_parameter_jacobian_type(model_types...) = Nonlinear()
coupling_time_gradient_type(model_types...) = Nonlinear()

isempty(::MatrixType) = false
isempty(::Empty) = true

iszero(::MatrixType) = false
iszero(::Zeros) = true

isidentity(::MatrixType) = false
isidentity(::Identity) = true

isinvariant(::MatrixType) = false
isinvariant(::Empty) = true
isinvariant(::Zeros) = true
isinvariant(::Identity) = true

isconstant(::MatrixType) = false
isconstant(::Empty) = true
isconstant(::Zeros) = true
isconstant(::Identity) = true
isconstant(::Invariant) = true
isconstant(::Constant) = true

islinear(::MatrixType) = false
islinear(::Empty) = true
islinear(::Zeros) = true
islinear(::Identity) = true
islinear(::Invariant) = true
islinear(::Constant) = true
islinear(::Linear) = true
