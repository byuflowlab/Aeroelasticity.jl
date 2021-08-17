"""
    LiftingLineFlaps{NF,NG,TF,TG} <: AbstractModel

Lifting line control surface model with `NF` cross sections and `NG` independent
control inputs, constructed by using the control surface models in `TF` and the
gains for each control input in `TG`.  State variables, inputs, and parameters
correspond to the state variables and parameters of each of the cross sections
concatenated.
"""
struct LiftingLineFlaps{NF,NG,TF<:NTuple{NF,Any},TG<:NTuple{NG,Any}} <: AbstractModel
    models::TF
    gains::TG
end

# --- Constructors --- #

"""
    LiftingLineFlaps(models::NTuple, gains::NTuple)

Construct a lifting line control surface model given the control surface models for each
section and the gains associated with each control input.
"""
LiftingLineFlaps(models, gains)

"""
    LiftingLineFlaps{NF}(model, gains)

Construct a lifting line control surface model using `NF` instances of `model`
and the gains in `gains` associated with each control input.
"""
function LiftingLineFlaps{NF}(model, gains) where NF
    models = ntuple(i->model, NF)
    return LiftingLineFlaps(models, gains)
end

# --- Traits --- #

function number_of_states(model::LiftingLineFlaps)
    return sum(number_of_states.(model.models))
end

function number_of_inputs(model::LiftingLineFlaps)
    return sum(number_of_inputs.(model.models))
end

number_of_parameters(model::LiftingLineFlaps) = sum(number_of_parameters.(model.models))

inplaceness(::Type{<:LiftingLineFlaps}) = InPlace()

function mass_matrix_type(::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NF,NG,TF,TG}
    model_types = (TF.parameters...,)
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

function state_jacobian_type(::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NF,NG,TF,TG}
    model_types = (TF.parameters...,)
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

function input_jacobian_type(::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NF,NG,TF,TG}
    model_types = (TF.parameters...,)
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


# --- Convenience Methods --- #

function set_states!(x, model::LiftingLineFlaps; section_states)

    section_models = model.models

    section_indices = state_indices(section_models)

    vxs = view.(Ref(x), section_indices)

    bfn! = (x, model, states) -> set_states!(x, model; states...)

    bfn!.(vxs, section_models, section_states)

    return x
end

function set_inputs!(y, model::LiftingLineFlaps; section_inputs)

    section_models = model.models

    section_indices = input_indices(section_models)

    vys = view.(Ref(y), section_indices)

    bfn! = (y, model, inputs) -> set_inputs!(y, model; inputs...)

    bfn!.(vys, section_models, section_inputs)

    return y
end

function set_parameters!(p, model::LiftingLineFlaps; section_parameters)

    section_models = model.models

    section_indices = parameter_indices(section_models)

    vps = view.(Ref(p), section_indices)

    bfn! = (p, model, parameters) -> set_parameters!(p, model; parameters...)

    bfn!.(vps, section_models, section_parameters)

    return p
end

function separate_states(model::LiftingLineFlaps, x)

    section_models = model.models

    section_indices = state_indices(section_models)

    vxs = view.(Ref(x), section_indices)

    return (section_states = separate_states.(section_models, vxs),)
end

function separate_inputs(model::LiftingLineFlaps, y)

    section_models = model.models

    section_indices = input_indices(section_models)

    vys = view.(Ref(y), section_indices)

    return (section_inputs = separate_inputs.(section_models, vys),)
end

function separate_parameters(model::LiftingLineFlaps, p)

    section_models = model.models

    section_indices = parameter_indices(section_models)

    vps = view.(Ref(p), section_indices)

    return (section_parameters = separate_parameters.(section_models, vps),)
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
