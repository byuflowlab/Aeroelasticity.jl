"""
    LiftingLineFlaps{NS,NF,T} <: AbstractModel

Lifting line control surface model with `NS` cross sections and `NF`
independent flaps, constructed by using the control surface models in `T`.
State variables and parameters correspond to the state variables and parameters
of each of the cross sections concatenated.  Inputs correspond to the state
variables and parameters of each of the cross sections concatenated, followed by
the deflection angle for each control surface.
"""
struct LiftingLineFlaps{NS,NF,T} <: AbstractModel
    models::T
    flaps::NTuple{NF, Vector{Int}}
end

"""
    LiftingLineSectionControl <: AbstractModel

Lifting line flap section model with state variables ``\\delta_1, \\delta_2,
\\dots, \\delta_N``, zero inputs, and zero parameters.  Two-dimensional control
surface models may be extended to three dimensional models by coupling with this
model.  Note that this model has no rate equations of its own since its state
variables are defined as functions of the 3D system's control variables.
"""
struct LiftingLineSectionControl <: AbstractModel end

number_of_states(::Type{LiftingLineSectionControl{NS,NF,T}}) where {NS,NF,T} = NF
number_of_inputs(::Type{<:LiftingLineSectionControl}) = 0
number_of_parameters(::Type{<:LiftingLineSectionControl}) = 0
inplaceness(::Type{LiftingLineSectionControl}) = OutOfPlace()

# --- Constructors --- #

"""
    LiftingLineFlaps(models, flaps)

Construct a lifting line control surface model given a tuple of 2D control surface
models corresponding to each section and the control deflection gains for each
control surface.
"""
function LiftingLineFlaps(models::T, flaps::TF) where {T<:NTuple{NS,Any},TF<:NTuple{NF,Any}} where {NS,NF}
    return LiftingLineFlaps{NS,NF,T}(models, flaps)
end

function LiftingLineFlaps{NS}(models::T, flaps::TF) where {T<:NTuple{NS,Any},TF<:NTuple{NF,Any}} where {NS,NF}
    return LiftingLineFlaps{NS,NF,T}(models, flaps)
end

function LiftingLineFlaps{NS,NF}(models::T, flaps::TF) where {T<:NTuple{NS,Any},TF<:NTuple{NF,Any}} where {NS,NF}
    return LiftingLineFlaps{NS,NF,T}(models, flaps)
end

"""
    LiftingLineFlaps{NS}(model, flaps)

Construct a lifting line control surface model using `NS` instances of `model`
and the control deflection gains for each control surface.
"""
function LiftingLineFlaps{NS}(model, flaps) where {N,T}
    models = ntuple(i->model, N)
    return LiftingLineFlaps(models, flaps)
end

# --- Traits --- #

function number_of_states(model::LiftingLineFlaps)
    return sum(number_of_states.(model.models))
end

function number_of_inputs(model::LiftingLineFlaps{NS,NF,T}) where {NS,NF,T}
    return sum(number_of_inputs.(model.models)) + NF
end

number_of_parameters(model::LiftingLineFlaps) = sum(number_of_parameters.(model.models))

inplaceness(::Type{<:LiftingLineFlaps}) where  = InPlace()

function mass_matrix_type(::Type{LiftingLineFlaps{NS,NF,T}}) where {NS,NF,T}
    model_types = (T.parameters...,)
    if all(isempty.(mass_matrix_type.(model_types)))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types)))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types)))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLineFlaps{NS,NF,T}}) where {NS,NF,T}
    model_types = (T.parameters...,)
    if all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function input_jacobian_type(::Type{LiftingLineFlaps{NS,NF,T}}) where {NS,NF,T}
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

# --- Methods --- #

function get_rates!(du, model::LiftingLineFlaps, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dus = view.(Ref(du), iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_rates!.(dus, models, us, ys, ps, t)

    return du
end

function get_mass_matrix!(M, model::LiftingLineFlaps)

    M .= 0

    models = model.models

    iu = state_indices(model.models)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models)

    return M
end

function get_mass_matrix!(M, model::LiftingLineFlaps, u, y, p, t)

    M .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Ms = view.(Ref(M), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_mass_matrix!.(Ms, models, us, ys, ps, t)

    return M
end

# --- Performance Overloads --- #

function get_state_jacobian!(J, model::LiftingLineFlaps)

    J .= 0

    models = model.models

    iu = state_indices(model.models)

    Js = view.(Ref(J), iu, iu)

    get_state_jacobian!.(Js, models)

    return J
end

function get_state_jacobian!(J, model::LiftingLineFlaps, u, y, p, t)

    J .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_state_jacobian!.(Js, models, us, ys, ps, t)

    return J
end

function get_input_jacobian(model::LiftingLineFlaps)

    f! = (y, x) -> input_jacobian_product!(y, x, model)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

function get_input_jacobian(model::LiftingLineFlaps, λ, d, p, t)

    f! = (y, x) -> input_jacobian_product!(y, x, model, λ, d, p, t)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

# --- Unit Testing Methods --- #

function get_lhs(model::LiftingLineFlaps, du, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dus = view.(Ref(du), iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    return vcat(get_lhs.(models, dus, us, ys, ps, t)...)
end

# --- Internal Methods --- #

function input_jacobian_product!(y, x, model::LiftingLineFlaps)

    models = model.models

    iy = state_indices(models)
    ix = input_indices(models)

    Jyi = get_input_jacobian.(models)

    yi = view.(Ref(y), iy)
    xi = view.(Ref(x), ix)

    mul!.(yi, Jyi, xi)

    return y
end

function input_jacobian_product!(y, x, model::LiftingLineFlaps, λ, d, p, t)

    models = model.models

    iu = state_indices(models)
    id = input_indices(models)
    ip = parameter_indices(models)

    xi = view.(Ref(x), id)
    yi = view.(Ref(y), iu)

    λi = view.(Ref(λ), iu)
    di = view.(Ref(d), id)
    pi = view.(Ref(p), ip)

    Jyi = get_input_jacobian.(models, λi, di, pi, t)

    mul!.(yi, Jyi, xi)

    return y
end
