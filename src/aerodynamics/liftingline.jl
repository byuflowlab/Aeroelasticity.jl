"""
    LiftingLine{NA,TA} <: AbstractModel

Lifting line model with `NA` cross sections, using the aerodynamic models in `TA`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated
"""
struct LiftingLine{NA,TA<:NTuple{NA,Any}} <: AbstractModel
    models::TA
end

# --- Constructors --- #

"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of 2D aerodynamic models.
"""
LiftingLine(models)

LiftingLine{NA}(models::TA) where TA<:NTuple{NA,Any} where NA = LiftingLine{NA,TA}(models)

"""
    LiftingLine{NA}(model)

Construct a lifting line aerodynamic model using `NA` instances of `model`.
"""
LiftingLine{NA}(model) where NA = LiftingLine(ntuple(i -> model, NA))

# --- Traits --- #

number_of_states(model::LiftingLine) = sum(number_of_states.(model.models))
number_of_inputs(model::LiftingLine) = sum(number_of_inputs.(model.models))
number_of_parameters(model::LiftingLine) = sum(number_of_parameters.(model.models))

inplaceness(::Type{<:LiftingLine}) = InPlace()

function rate_jacobian_type(::Type{LiftingLine{NA,TA}}) where {NA,TA}
    model_types = (TA.parameters...,)
    if all(isempty.(rate_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(rate_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(rate_jacobian_type.(model_types)))
        return Identity()
    elseif all(isinvariant.(rate_jacobian_type.(model_types)))
        return Invariant()
    elseif all(isconstant.(rate_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(rate_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{NA,TA}}) where {NA,TA}
    model_types = (TA.parameters...,)
    if all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif all(isinvariant.(state_jacobian_type.(model_types)))
        return Invariant()
    elseif all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function input_jacobian_type(::Type{LiftingLine{NA,TA}}) where {NA,TA}
    model_types = (TA.parameters...,)
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

function parameter_jacobian_type(::Type{LiftingLine{NA,TA}}) where {NA,TA}
    model_types = (TA.parameters...,)
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

function time_gradient_type(::Type{LiftingLine{NA,TA}}) where {NA,TA}
    model_types = (TA.parameters...,)
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

# --- Methods --- #

function get_residual!(resid, model::LiftingLine, dx, x, y, p, t)

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    resids = view.(Ref(resid), ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_residual!.(resids, models, dxs, xs, ys, ps, t)

    return resid
end

# --- Performance 0verloads --- #

function get_rate_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    ix = state_indices(model.models)

    Js = view.(Ref(J), ix, ix)

    get_rate_jacobian!.(Js, models)

    return J
end

function get_rate_jacobian!(J, model::LiftingLine, p)

    J .= 0

    models = model.models

    ix = state_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)

    get_rate_jacobian!.(Js, models, ps)

    return J
end

function get_rate_jacobian!(J, model::LiftingLine, dx, x, y, p, t)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_rate_jacobian!.(Js, models, dxs, xs, ys, ps, t)

    return J
end

function get_state_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    ix = state_indices(model.models)

    Js = view.(Ref(J), ix, ix)

    get_state_jacobian!.(Js, models)

    return J
end

function get_state_jacobian!(J, model::LiftingLine, p)

    J .= 0

    models = model.models

    ix = state_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)

    get_state_jacobian!.(Js, models, ps)

    return J
end

function get_state_jacobian!(J, model::LiftingLine, dx, x, y, p, t)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_state_jacobian!.(Js, models, dxs, xs, ys, ps, t)

    return J
end

function get_input_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)

    Js = view.(Ref(J), ix, iy)

    get_input_jacobian!.(Js, models)

    return J
end

function get_input_jacobian!(J, model::LiftingLine, p)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, iy)
    ps = view.(Ref(p), ip)

    get_input_jacobian!.(Js, models, ps)

    return J
end

function get_input_jacobian!(J, model::LiftingLine, dx, x, y, p, t)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, iy)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_input_jacobian!.(Js, models, dxs, xs, ys, ps, t)

    return J
end

function get_parameter_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    ix = state_indices(model.models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ip)

    get_parameter_jacobian!.(Js, models)

    return J
end

function get_parameter_jacobian!(J, model::LiftingLine, p)

    J .= 0

    models = model.models

    ix = state_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ip)
    ps = view.(Ref(p), ip)

    get_parameter_jacobian!.(Js, models, ps)

    return J
end

function get_parameter_jacobian!(J, model::LiftingLine, dx, x, y, p, t)

    J .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), ix, ip)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_parameter_jacobian!.(Js, models, dxs, xs, ys, ps, t)

    return J
end

function get_time_gradient!(dT, model::LiftingLine)

    dT .= 0

    models = model.models

    ix = state_indices(model.models)

    dTs = view.(Ref(dT), ix)

    get_time_gradient!.(dTs, models)

    return dT
end

function get_time_gradient!(dT, model::LiftingLine, p)

    dT .= 0

    models = model.models

    ix = state_indices(models)
    ip = parameter_indices(models)

    dTs = view.(Ref(dT), ix)
    ps = view.(Ref(p), ip)

    get_time_gradient!.(dTs, models, ps)

    return dT
end

function get_time_gradient!(dT, model::LiftingLine, dx, x, y, p, t)

    dT .= 0

    models = model.models

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dTs = view.(Ref(dT), ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_time_gradient!.(dTs, models, dxs, xs, ys, ps, t)

    return dT
end

# --- Convenience Methods --- #

function set_states!(x, model::LiftingLine; section_states)

    section_models = model.models

    section_indices = state_indices(section_models)

    vxs = view.(Ref(x), section_indices)

    bfn! = (x, model, states) -> set_states!(x, model; states...)

    bfn!.(vxs, section_models, section_states)

    return x
end

function set_inputs!(y, model::LiftingLine; section_inputs)

    section_models = model.models

    section_indices = input_indices(section_models)

    vys = view.(Ref(y), section_indices)

    bfn! = (y, model, inputs) -> set_inputs!(y, model; inputs...)

    bfn!.(vys, section_models, section_inputs)

    return y
end

function set_parameters!(p, model::LiftingLine; section_parameters)

    section_models = model.models

    section_indices = parameter_indices(section_models)

    vps = view.(Ref(p), section_indices)

    bfn! = (p, model, parameters) -> set_parameters!(p, model; parameters...)

    bfn!.(vps, section_models, section_parameters)

    return p
end

function separate_states(model::LiftingLine, x)

    section_models = model.models

    section_indices = state_indices(section_models)

    vxs = view.(Ref(x), section_indices)

    return (section_states = separate_states.(section_models, vxs),)
end

function separate_inputs(model::LiftingLine, y)

    section_models = model.models

    section_indices = input_indices(section_models)

    vys = view.(Ref(y), section_indices)

    return (section_inputs = separate_inputs.(section_models, vys),)
end

function separate_parameters(model::LiftingLine, p)

    section_models = model.models

    section_indices = parameter_indices(section_models)

    vps = view.(Ref(p), section_indices)

    return (section_parameters = separate_parameters.(section_models, vps),)
end

# --- Internal Methods --- #
