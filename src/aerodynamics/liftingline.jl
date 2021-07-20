"""
    LiftingLine{N,T} <: AbstractModel

Lifting line model with `N` cross sections, using the aerodynamic models in `T`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated
"""
struct LiftingLine{N,T} <: AbstractModel
    models::T
end

"""
    LiftingLineSection <: AbstractModel

Lifting line section model with state variables ``q = \\begin{bmatrix} v_x & v_y &
v_z & \\omega_x & \\omega_y & \\omega_z \\end{bmatrix}^T``, inputs
``r = \\begin{bmatrix} L' & Y' & D' & Mx, My, Mz \\end{bmatrix}^T``, and no
parameters.  Two-dimensional aerodynamic models may be extended to
three-dimensional models by coupling with this model.  Note that
this model has no rate equations of its own since its state variables are
defined as functions of the 3D structural model's state variables.
"""
struct LiftingLineSection <: AbstractModel end

number_of_states(::Type{<:LiftingLineSection}) = 6
number_of_inputs(::Type{<:LiftingLineSection}) = 6
number_of_parameters(::Type{<:LiftingLineSection}) = 0
inplaceness(::Type{LiftingLineSection}) = OutOfPlace()

# --- Constructors --- #

"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of aerodynamic models.
"""
LiftingLine(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

LiftingLine{N}(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

"""
    LiftingLine{N}(model)

Construct a lifting line aerodynamic model using `N` instances of `model`.
"""
function LiftingLine{N}(model) where {N,T}
    models = ntuple(i->model, N)
    return LiftingLine(models)
end

# --- Traits --- #

number_of_states(model::LiftingLine) = sum(number_of_states.(model.models))
number_of_inputs(model::LiftingLine) = sum(number_of_inputs.(model.models))
number_of_parameters(model::LiftingLine) = sum(number_of_parameters.(model.models))
inplaceness(::Type{LiftingLine{N,T}}) where {N,T} = InPlace()

function mass_matrix_type(::Type{LiftingLine{N,T}}) where {N,T}
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

function state_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
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

function input_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
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

function get_rates!(du, model::LiftingLine, u, y, p, t)

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

function get_mass_matrix!(M, model::LiftingLine)

    M .= 0

    models = model.models

    iu = state_indices(model.models)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models)

    return M
end

function get_mass_matrix!(M, model::LiftingLine, u, y, p, t)

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

function get_state_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    iu = state_indices(model.models)

    Js = view.(Ref(J), iu, iu)

    get_state_jacobian!.(Js, models)

    return J
end

function get_state_jacobian!(J, model::LiftingLine, u, y, p, t)

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

function get_input_jacobian(model::LiftingLine)

    f! = (y, x) -> input_jacobian_product!(y, x, model)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

function get_input_jacobian(model::LiftingLine, λ, d, p, t)

    f! = (y, x) -> input_jacobian_product!(y, x, model, λ, d, p, t)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

# --- Unit Testing Methods --- #

function get_lhs(model::LiftingLine, du, u, y, p, t)

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

function input_jacobian_product!(y, x, model::LiftingLine)

    models = model.models

    iy = state_indices(models)
    ix = input_indices(models)

    Jyi = get_input_jacobian.(models)

    yi = view.(Ref(y), iy)
    xi = view.(Ref(x), ix)

    mul!.(yi, Jyi, xi)

    return y
end

function input_jacobian_product!(y, x, model::LiftingLine, λ, d, p, t)

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
