"""
    LiftingLine{N,T} <: AbstractModel

Lifting line model with `N` cross sections, using the aerodynamic models in `T`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated.
"""
struct LiftingLine{N,T} <: AbstractModel
    models::T
end

# --- Constructors --- #
"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of aerodynamic models.
"""
LiftingLine(models::NTuple{N,T}) where {N,T} = LiftingLine{N}(models)

# --- Traits --- #
number_of_states(model::LiftingLine) = number_of_states(model.models)
number_of_inputs(model::LiftingLine) = number_of_inputs(model.models)
number_of_parameters(model::LiftingLine) = number_of_parameters(model.models)
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
    else
        return Varying()
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
    else
        return Varying()
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
    else
        return Varying()
    end
end

function input_dependence_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(_linear_input_dependence.(model_types))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- Methods --- #

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

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models, us, ys, ps)

    return M
end

function get_rates!(du, model::LiftingLine, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    dus = view.(Ref(du), iu)

    get_rates!.(dus, models, us, ys, ps)

    return du
end

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

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Js = view.(Ref(J), iu, iu)

    get_state_jacobian!.(Js, models, us, ys, ps)

    return J
end

function get_input_jacobian!(Jy, model::LiftingLine)

    models = model.models

    iu = state_indices(model.models)
    iy = input_indices(model.models)

    Jys = view.(Ref(Jy), iu, iy)

    get_state_jacobian!.(Jys, models)

    return Jy
end

function get_input_jacobian!(Jy, model::LiftingLine, u, y, p, t)

    Jy .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Jys = view.(Ref(Jy), iu, iy)

    get_state_jacobian!.(Jys, models, us, ys, ps)

    return Jy
end

# TODO: Add parameter jacobian
