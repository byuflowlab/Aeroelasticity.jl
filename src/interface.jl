const number_of_rates = number_of_states

"""
    number_of_states(model)

Return the total number of state variables corresponding to a model.
"""
number_of_states

number_of_states(model::Model) = val(model.nx)

number_of_states(coupling::Coupling) = val(coupling.nx)

number_of_states(model::CoupledModel) = number_of_states(model.coupling)

"""
    number_of_inputs(model)

Return the total number of inputs corresponding to a model.
"""
number_of_inputs

number_of_inputs(model::Model) = val(model.ny)

number_of_inputs(coupling::Coupling) = val(coupling.ny)

number_of_inputs(model::CoupledModel) = number_of_inputs(model.coupling)

"""
    number_of_parameters(model)

Return the total number of parameters corresponding to a model.
"""
number_of_parameters

number_of_parameters(model::Model) = val(model.np)

number_of_parameters(coupling::Coupling) = val(coupling.np)

number_of_parameters(model::CoupledModel) = number_of_parameters(model.coupling)

"""
    number_of_additional_parameters(model)

Return the total number of additional parameters corresponding to a coupled model.
"""
number_of_additional_parameters

number_of_additional_parameters(coupling::Coupling) = val(coupling.npc)

number_of_additional_parameters(model::CoupledModel) = number_of_additional_parameters(model.coupling)

const rate_indices = state_indices

"""
    state_indices(model::CoupledModel)

Return the state variable indices for each submodel in `model`
"""
state_indices(model::CoupledModel) = state_indices(model.models)

function state_indices(models)
    nx = number_of_states.(models)
    ix2 = cumsum(nx)
    ix1 = ix2 .- nx .+ 1
    return UnitRange.(ix1, ix2)
end

"""
    input_indices(model::CoupledModel)

Return the input variable indices for each submodel in `model`
"""
input_indices(model::CoupledModel) = input_indices(model.models)

function input_indices(models)
    ny = number_of_inputs.(models)
    iy2 = cumsum(ny)
    iy1 = iy2 .- ny .+ 1
    return UnitRange.(iy1, iy2)
end

"""
    parameter_indices(model::CoupledModel)

Return the parameter indices for each submodel in `model`
"""
parameter_indices(model::CoupledModel) = parameter_indices(model.models)

function parameter_indices(models)
    np = number_of_parameters.(models)
    ip2 = cumsum(np)
    ip1 = ip2 .- np .+ 1
    return UnitRange.(ip1, ip2)
end

"""
    additional_parameter_indices(model::CoupledModel)

Return the indices corresponding to additional parameters introduced by the model coupling
"""
function additional_parameter_indices(model::CoupledModel)
    submodels = model.models
    ip1 = sum(number_of_parameters.(submodels)) + 1
    ip2 = number_of_parameters(model)
    return UnitRange(ip1, ip2)
end

const get_rates = get_states

"""
    get_states(model::Model; kwargs...)

Return the state vector for `model` using the values in `kwargs`.
"""
function get_states(model::Model; kwargs...)
    x = zeros(number_of_states(model))
    set_states!(x, model; kwargs...)
    return x
end

"""
    get_inputs(model::Model; kwargs...)

Return the input vector for `model` using the values in `kwargs`.
"""
function get_inputs(model::Model; kwargs...)
    y = zeros(number_of_inputs(model))
    set_inputs!(y, model; kwargs...)
    return y
end

"""
    get_parameters(model::Model; kwargs...)

Return the parameter vector for `model` using the values in `kwargs`.
"""
function get_parameters(model::Model; kwargs...)
    p = zeros(number_of_parameters(model))
    set_parameters!(p, model; kwargs...)
    return p
end

function get_parameters(coupling::Coupling; kwargs...)
    p = zeros(number_of_additional_parameters(coupling))
    set_parameters!(p, coupling; kwargs...)
    return p
end

"""
    get_additional_parameters(model::CoupledModel; kwargs...)

Return the additional parameter vector for `model` using the values in `kwargs`.
"""
function get_additional_parameters(model::CoupledModel; kwargs...)
    return get_parameters(model.coupling; kwargs...)
end

const set_rates! = set_states!

"""
    set_states!(x, model::Model; kwargs...)

In-place version of [`get_states`](@ref)
"""
set_states!(x, model::Model; kwargs...) = model.setstate(x; kwargs...)

"""
    set_inputs!(y, model::Model; kwargs...)

In-place version of [`get_inputs`](@ref)
"""
set_inputs!(y, model::Model; kwargs...) = model.setinput(y; kwargs...)

"""
    set_parameters!(p, model::Model; kwargs...)

In-place version of [`get_parameters`](@ref)
"""
set_parameters!(p, model::Model; kwargs...) = model.setparam(p; kwargs...)

set_parameters!(p, coupling::Coupling; kwargs...) = coupling.setparam(p; kwargs...)

"""
    set_states!(x, model::CoupledModel, i; kwargs...)

Set the elements of the state vector `x` corresponding to the `i`th submodel of `model` 
to the values in `kwargs`
"""
function set_states!(x, model::CoupledModel, i; kwargs...)

    return set_states!(view(x, state_indices(model)[i]), model.models[i]; kwargs...)
end

"""
    set_inputs!(y, model::CoupledModel, i; kwargs...)

Set the elements of the input vector `y` corresponding to the `i`th submodel of `model` 
to the values in `kwargs`
"""
function set_inputs!(y, model::CoupledModel, i; kwargs...)

    return set_inputs!(view(y, input_indices(model)[i]), model.models[i]; kwargs...)
end

"""
    set_parameters!(p, model::CoupledModel, i; kwargs...)

Set the elements of the parameter vector `p` corresponding to the `i`th submodel of `model` 
to the values in `kwargs`
"""
function set_parameters!(p, model::CoupledModel, i; kwargs...)

    return set_parameters!(view(p, parameter_indices(model)[i]), model.models[i]; kwargs...)
end

"""
    set_additional_parameters!(p, model::CoupledModel; kwargs...)

Set the elements of the parameter vector `p` corresponding to additional parameters 
introduced by the model coupling to the values specified in `kwargs`
"""
function set_additional_parameters!(p, model::CoupledModel; kwargs...)

    return set_parameters!(view(p, additional_parameter_indices(model)), model.coupling; kwargs...)
end

const separate_rates = separate_states

"""
    separate_states(model, x)

Separate the state vector entries in `x`.
"""
separate_states

separate_states(model::Model, x) = model.sepstate(x)

function separate_states(model::CoupledModel, x)

    xs = view.(Ref(x), state_indices(model))

    return separate_states.(model.models, xs)
end

"""
    separate_inputs(model, y)

Separate the input vector entries in `y`
"""
separate_inputs

separate_inputs(model::Model, y) = model.sepinput(y)

function separate_inputs(model::CoupledModel, y)

    ys = view.(Ref(y), input_indices(model))

    return separate_inputs.(model.models, ys)
end

"""
    separate_parameters(model, p)

Separate the parameter vector entries in `p`
"""
separate_parameters

separate_parameters(model::Model, p) = model.sepparam(p)

separate_parameters(model::Coupling, p) = model.sepparam(p)

function separate_parameters(model::CoupledModel, p) where N

    pmodels = view.(Ref(p), parameter_indices(model))
    pmodels_sep = separate_parameters.(model.models, pmodels)

    padd = view(p, additional_parameter_indices(model))
    padd_sep = separate_parameters(model.coupling, padd)

    return (pmodels_sep..., padd_sep)
end

"""
    get_residual(model, dx, x, y, p, t)

Calculate the residual for the specified model.
"""
function get_residual(model, dx, x, y, p, t)

    return _get_residual(inplaceness(model), model, dx, x, y, p, t)
end

# dispatch to an in-place function
function _get_residual(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    resid = similar(dx, TF)
    get_residual!(resid, model, dx, x, y, p, t)
    return resid
end

# calculate residual for a standalone model
_get_residual(::OutOfPlace, model::Model, dx, x, y, p, t) = model.f(dx, x, y, p, t)

# calculate residual for a coupled model
function _get_residual(::OutOfPlace, model::CoupledModel, dx, x, y, p, t)

    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    return vcat(get_residual.(model.models, dxs, xs, ys, ps, t)...)
end

"""
    get_residual!(resid, model, dx, x, y, p, t)

In-place version of [`get_residual`](@ref)
"""
function get_residual!(resid, model, dx, x, y, p, t)
    return _get_residual!(resid, inplaceness(model), model, dx, x, y, p, t)
end

# dispatch to an out-of-place function
function _get_residual!(resid, ::OutOfPlace, model, dx, x, y, p, t)
    return resid .= get_residual(model, dx, x, y, p, t)
end

# calculate residual for a standalone model
function _get_residual!(resid, ::InPlace, model::Model, dx, x, y, p, t)
   return model.f(resid, dx, x, y, p, t)
end

# calculate residual for a coupled model
function _get_residual!(resid, ::InPlace, model::CoupledModel, dx, x, y, p, t)

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    resids = view(Ref(resid), ix[i])
    dxs = view.(Ref(dx), ix[i])
    xs = view.(Ref(xs), ix[i])
    ys = view.(Ref(ys), iy[i])
    ps = view.(Ref(ps), ip[i])

    get_residual!.(resids, models, dxs, xs, ys, ps, t)

    return resid
end

"""
    get_rate_jacobian(model)
    get_rate_jacobian(model, p)
    get_rate_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the residual function for `model` with respect to the
state rates
"""
function get_rate_jacobian(model, args...)
    return _get_rate_jacobian(inplaceness(model), model, args...)
end

# dispatch to an in-place function
function _get_rate_jacobian(::InPlace, model)
    Nx = number_of_states(model)
    J = Matrix{Float64}(undef, Nx, Nx)
    return get_rate_jacobian!(J, model)
end

# dispatch to an in-place function
function _get_rate_jacobian(::InPlace, model, p)
    TF = promote_type(eltype(p))
    Nx = number_of_states(model)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_rate_jacobian!(J, model, p)
end

# dispatch to an in-place function
function _get_rate_jacobian(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(model)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_rate_jacobian!(J, model, dx, x, y, p, t)
end

# calculate rate jacobian for a standalone model
function _get_rate_jacobian(::OutOfPlace, model::Model, args...)
    return jacobian(model.ratejac, model.nx, model.nx, args...)
end

# calculate rate jacobian for a coupled model
function _get_rate_jacobian(::OutOfPlace, model::CoupledModel)

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(model)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(model)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(model)

    # calculate jacobian (chain rule)
    return r_dx + r_y*y_dx
end

# calculate rate jacobian for a coupled model
function _get_rate_jacobian(::OutOfPlace, model::CoupledModel, p)

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(model, p)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(model, p)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(model, p)

    # calculate jacobian (chain rule)
    return r_dx + r_y*y_dx
end

# calculate rate jacobian for a coupled model
function _get_rate_jacobian(::OutOfPlace, model::CoupledModel, dx, x, y, p, t)

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(model, dx, x, y, p, t)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(model, dx, x, y, p, t)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(model, dx, x, p, t)

    return r_dx + r_y*y_dx
end

function get_rate_partial(model::CoupledModel)

    # initialize diagonal blocks
    Jii = get_rate_jacobian.(model.models)

    # assemble into block diagonal matrix
    return block_diagonal(Jii)
end

function get_rate_partial(model::CoupledModel, p) where N
    
    # get parameters for each submodel
    ip = parameter_indices(model)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_rate_jacobian.(model.models, ps)

    # assemble into block diagonal matrix
    return block_diagonal(Jii)
end

function get_rate_partial(model::CoupledModel, dx, x, y, p, t) where N

    # get variable indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get views into variables for each submodel
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_rate_jacobian.(model.models, dxs, xs, ys, ps, t)

    # assemble into block diagonal matrix
    return block_diagonal(Jii)
end

"""
    get_rate_jacobian!(J, model)
    get_rate_jacobian!(J, model, p)
    get_rate_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_rate_jacobian`](@ref).
"""
function get_rate_jacobian!(J, model, args...; kwargs...)
    return _get_rate_jacobian!(J, inplaceness(model), model, args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_rate_jacobian!(J, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_rate_jacobian(model, args...; kwargs...)
end

# calculate rate jacobian for a standalone model
function _get_rate_jacobian!(J, ::InPlace, model::Model, args...)
    return jacobian!(J, model.ratejac, args...)
end

# calculate rate jacobian for a coupled model
function _get_rate_jacobian!(J, ::InPlace, model::CoupledModel;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dyddx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    )

    # get state and input indices
    ix = state_indices(model)
    iy = input_indices(model)

    # calculate input mass matrix
    get_coupling_rate_jacobian!(dyddx_cache, model)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, submodels[i])
        for j = 1:length(submodels)
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, submodels[i])
                mul!(ri_dxj, ri_yi, yi_dxj, 1, 1)
            else
                mul!(ri_dxj, ri_yi, yi_dxj)
            end
        end
    end

    return J
end

# calculate rate jacobian for a coupled model
function _get_rate_jacobian!(J, ::InPlace, model::CoupledModel, p;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dyddx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    )

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # calculate input rate jacobian
    get_coupling_rate_jacobian!(dyddx_cache, model, p)

    # calculate rate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        get_input_jacobian!(ri_yi, submodels[i], pi)
        for j = 1:length(submodels)
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, submodels[i], pi)
                mul!(ri_dxj, ri_yi, yi_dxj, 1, 1)
            else
                mul!(ri_dxj, ri_yi, yi_dxj)
            end
        end
    end

    return J
end

# calculate rate jacobian for a combination of models
function _get_rate_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace, model::CoupledModel, 
    dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dyddx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    ) where N

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # calculate input rate jacobian
    get_coupling_rate_jacobian!(dyddx_cache, model, dx, x, p, t)

    # calculate rate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, submodels[i], dxi, xi, yi, pi, t)
        for j = 1:length(submodels)
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, submodels[i], dxi, xi, yi, pi, t)
                mul!(ri_dxj, ri_yi, yi_dxj, 1, 1)
            else
                mul!(ri_dxj, ri_yi, yi_dxj)
            end
        end
    end

    return J
end

"""
    get_state_jacobian(model)
    get_state_jacobian(model, p)
    get_state_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the residual function for `model` with respect to the
state variables
"""
function get_state_jacobian(model, args...)
    return _get_state_jacobian(inplaceness(model), model, args...)
end

# dispatch to an in-place function
function _get_state_jacobian(::InPlace, model)
    Nx = number_of_states(model)
    J = Matrix{Float64}(undef, Nx, Nx)
    return get_state_jacobian!(J, model)
end

# dispatch to an in-place function
function _get_state_jacobian(::InPlace, model, p)
    TF = eltype(p)
    Nx = number_of_states(model)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_state_jacobian!(J, model, p)
end

# dispatch to an in-place function
function _get_state_jacobian(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(model)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_state_jacobian!(J, model, dx, x, y, p, t)
end

# calculate state jacobian for a standalone model
function _get_state_jacobian(::OutOfPlace, model::Model, args...)
    return jacobian(model.statejac, model.nx, model.nx, args...)
end

# calculate state jacobian for a coupled model
function _get_state_jacobian(::OutOfPlace, model::CoupledModel) where N

    # initialize mass matrix
    r_x = get_state_partial(model)

    # calculate input jacobian
    r_y = get_input_jacobian(model)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(model)

    return r_x + r_y*y_x
end

# calculate state jacobian for a coupled model
function _get_state_jacobian(::OutOfPlace, model::CoupledModel, p) where N

    # initialize mass matrix
    r_x = get_state_partial(model, p)

    # calculate input jacobian
    r_y = get_input_jacobian(model, p)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(model, p)

    return r_x + r_y*y_x
end

# calculate state jacobian for a coupled model
function _get_state_jacobian(::OutOfPlace, model::CoupledModel, dx, x, y, p, t) where N

    # initialize mass matrix
    r_x = get_state_partial(model, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(model, dx, x, y, p, t)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(model, dx, x, p, t)

    return r_x + r_y*y_x
end

function get_state_partial(model::CoupledModel) where N

    submodels = model.models

    # initialize diagonal blocks
    Jii = get_state_jacobian.(submodels)

    # assemble block diagonal jacobian
    return block_diagonal(Jii)
end

function get_state_partial(model::CoupledModel, p) where N
    
    submodels = model.models

    ip = parameter_indices(model)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_state_jacobian.(submodels, ps)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

function get_state_partial(model::CoupledModel, dx, x, y, p, t) where N
    
    submodels = model.models

    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_state_jacobian.(submodels, dxs, xs, ys, ps, t)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

"""
    get_state_jacobian!(J, model)
    get_state_jacobian!(J, model, p)
    get_state_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_state_jacobian`](@ref)
"""
function get_state_jacobian!(J, model, args...; kwargs...)
    return _get_state_jacobian!(J, inplaceness(model), model, args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_state_jacobian!(J, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_state_jacobian(model, args...; kwargs...)
end

# calculate state jacobian for a standalone model
function _get_state_jacobian!(J, ::InPlace, model::Model, args...)
    return jacobian!(J, model.statejac, args...)
end

# calculate state jacobian for a coupled model
function _get_state_jacobian!(J, ::InPlace, model::CoupledModel;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    ) where N

    # get state and input indices
    ix = state_indices(model)
    iy = input_indices(model)

    # get input jacobian
    get_coupling_state_jacobian!(dydx_cache, model)

    # calculate jacobian
    for i = 1:length(model.models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, model.models[i])
        for j = 1:length(model.models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, model.models[i])
                mul!(ri_xj, ri_yi, yi_xj, 1, 1)
            else
                mul!(ri_xj, ri_yi, yi_xj)
            end
        end
    end

    return J
end

# calculate state jacobian for a coupled model
function _get_state_jacobian!(J, ::InPlace, model::CoupledModel, p;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    ) where N

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input state jacobian
    get_coupling_state_jacobian!(dydx_cache, model, p)

    # calculate jacobian
    for i = 1:length(model.models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, model.models[i], pi)
        for j = 1:length(model.models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, model.models[i], pi)
                mul!(ri_xj, ri_yi, yi_xj, 1, 1)
            else
                mul!(ri_xj, ri_yi, yi_xj)
            end
        end
    end

    return J
end

# calculate state jacobian for a coupled model
function _get_state_jacobian!(J, ::InPlace, model::CoupledModel, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydx_cache = similar(J, number_of_inputs(model), number_of_states(model))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input state jacobian
    get_coupling_state_jacobian!(dydx_cache, model, dx, x, p, t)

    # calculate jacobian
    for i = 1:length(model.models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, model.models[i], dxi, xi, yi, pi, t)
        for j = 1:length(model.models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, model.models[i], dxi, xi, yi, pi, t)
                mul!(ri_xj, ri_yi, yi_xj, 1, 1)
            else
                mul!(ri_xj, ri_yi, yi_xj)
            end
        end
    end

    return J
end

"""
    get_input_jacobian(model)
    get_input_jacobian(model, p)
    get_input_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the residual function for `model` with respect to the
input variables
"""
function get_input_jacobian(model, args...)
    return _get_input_jacobian(inplaceness(model), model, args...)
end

# dispatch to an in-place function
function _get_input_jacobian(::InPlace, model)
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    J = Matrix{Float64}(undef, Nx, Ny)
    return get_input_jacobian!(J, model)
end

# dispatch to an in-place function
function _get_input_jacobian(::InPlace, model, p)
    TF = eltype(p)
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    J = Matrix{TF}(undef, Nx, Ny)
    return get_input_jacobian!(J, model, p)
end

# dispatch to an in-place function
function _get_input_jacobian(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    J = Matrix{TF}(undef, Nx, Ny)
    return get_input_jacobian!(J, model, dx, x, y, p, t)
end

# calculate input jacobian for a standalone model
function _get_input_jacobian(::OutOfPlace, model::Model, args...)
    return jacobian(model.inputjac, model.nx, model.ny, args...)
end

# calculate input jacobian for a coupled model
function _get_input_jacobian(::OutOfPlace, model::CoupledModel) where N

    submodels = model.models

    # initialize diagonal blocks
    Jii = get_input_jacobian.(submodels)

    # assemble block diagonal jacobian
    return block_diagonal(Jii)
end

# calculate input jacobian for a coupled model
function _get_input_jacobian(::OutOfPlace, model::CoupledModel, p) where N

    submodels = model.models

    ip = parameter_indices(model)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_input_jacobian.(submodels, ps)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

# calculate input jacobian for a coupled model
function _get_input_jacobian(::OutOfPlace, model::CoupledModel, dx, x, y, p, t) where N

    submodels = model.models

    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_input_jacobian.(submodels, dxs, xs, ys, ps, t)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

"""
    get_input_jacobian!(J, model)
    get_input_jacobian!(J, model, p)
    get_input_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_input_jacobian`](@ref)
"""
function get_input_jacobian!(J, model, args...; kwargs...)
    return _get_input_jacobian!(J, inplaceness(model), model, args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_input_jacobian!(J, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_input_jacobian(model, args...; kwargs...)
end

# calculate input jacobian for a standalone model
function _get_input_jacobian!(J, ::InPlace, model::Model, args...)
    return jacobian!(J, model.inputjac, args...)
end

# calculate input jacobian for a coupled model
function _get_input_jacobian!(J, ::InPlace, model::CoupledModel)

    submodels = model.models

    # get state and input indices
    ix = state_indices(model)
    iy = input_indices(model)

    # calculate jacobian
    J .= 0
    for i = 1:length(submodels)
        Jii = view(J, ix[i], iy[j])
        get_state_jacobian!(Jii, submodels[i])
    end

    return J
end

# calculate input jacobian for a coupled model
function _get_input_jacobian!(J, ::Constant, ::InPlace, model::CoupledModel, p)

    submodels = model.models

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # calculate jacobian
    J .= 0
    for i = 1:length(submodels)
        Jii = view(J, ix[i], iy[j])
        pi = view(p, ip[i])
        get_input_jacobian!(Jii, submodels[i], pi)
    end

    return J
end

# calculate input jacobian for a coupled model
function _get_input_jacobian!(J, ::InPlace, model::CoupledModel, dx, x, y, p, t)

    submodels = model.models

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # calculate jacobian
    J .= 0
    for i = 1:length(submodels)
        Jii = view(J, ix[i], iy[j])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        get_input_jacobian!(Jii, submodels[i], dxi, xi, yi, pi, t)
    end

    return J
end

"""
    get_parameter_jacobian(model)
    get_parameter_jacobian(model, p)
    get_parameter_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the residual function for `model` with respect to the
parameters
"""
function get_parameter_jacobian(model, args...)
    return _get_parameter_jacobian(inplaceness(model), model, args...)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::InPlace, model)
    Nx = number_of_states(model)
    Np = number_of_parameters(model)
    J = Matrix{Float64}(undef, Nx, Np)
    return get_parameter_jacobian!(J, model)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::InPlace, model, p)
    TF = eltype(p)
    Nx = number_of_states(model)
    Np = number_of_parameters(model)
    J = Matrix{TF}(undef, Nx, Np)
    return get_parameter_jacobian!(J, model, p)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(model)
    Np = number_of_parameters(model)
    J = Matrix{TF}(undef, Nx, Np)
    return get_parameter_jacobian!(J, model, dx, x, y, p, t)
end

# calculate parameter jacobian for a standalone model
function _get_parameter_jacobian(::OutOfPlace, model::Model, args...)
    return jacobian(model.paramjac, model.nx, model.np, args...)
end

# calculate parameter jacobian for a coupled model
function _get_parameter_jacobian(::OutOfPlace, model::CoupledModel)

    # initialize mass matrix
    r_p = get_parameter_partial(model)

    # calculate input jacobian
    r_y = get_input_jacobian(model)

    # calculate input jacobian
    y_p = get_coupling_parameter_jacobian(model)

    return r_p + r_y*y_p
end

# calculate parameter jacobian for a coupled model
function _get_parameter_jacobian(::OutOfPlace, model::CoupledModel, p)

    # initialize parameter jacobian
    r_p = get_parameter_partial(model, p)

    # calculate input jacobian
    r_y = get_input_jacobian(model, p)

    # calculate coupling function parameter jacobian
    y_p = get_coupling_parameter_jacobian(model, p)

    return r_p + r_y*y_p
end

# calculate parameter jacobian for a coupled model
function _get_parameter_jacobian(::OutOfPlace, model::CoupledModel, dx, x, y, p, t)

    # initialize parameter jacobian
    r_p = get_parameter_partial(model, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(model, dx, x, y, p, t)

    # calculate coupling function parameter jacobian
    y_p = get_coupling_parameter_jacobian(model, dx, x, p, t)

    return r_p + r_y*y_p
end

function get_parameter_partial(model::CoupledModel)

    submodels = model.models

    # initialize diagonal blocks
    Jii = get_parameter_jacobian.(submodels)

    # assemble block diagonal jacobian
    return block_diagonal(Jii)
end

function get_parameter_partial(model::CoupledModel, p)
    
    submodels = model.models

    ip = parameter_indices(model)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_parameter_jacobian.(submodels, ps)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

function get_parameter_partial(model::CoupledModel, dx, x, y, p, t)
    
    submodels = model.models

    ix = parameter_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    # initialize diagonal blocks
    Jii = get_parameter_jacobian.(submodels, dxs, xs, ys, ps, t)

    # return block diagonal matrix
    return block_diagonal(Jii)
end

"""
    get_parameter_jacobian!(J, model)
    get_parameter_jacobian!(J, model, p)
    get_parameter_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_parameter_jacobian`](@ref)
"""
function get_parameter_jacobian!(J, model, args...; kwargs...)
    return _get_parameter_jacobian!(J, inplaceness(model), model, args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_parameter_jacobian!(J, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_parameter_jacobian(model, args...; kwargs...)
end

# calculate parameter jacobian for a standalone model
function _get_parameter_jacobian!(J, ::InPlace, model::Model, args...)
    return jacobian!(J, model.paramjac, args...)
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian!(J, ::InPlace, model::CoupledModel;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydp_cache = similar(J, number_of_inputs(model), number_of_parameters(model)))

    # get state and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input parameter jacobian
    get_coupling_parameter_jacobian!(dydp_cache, model)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, submodels[i])
        for j = 1:length(submodels)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, submodels[i])
                mul!(ri_pj, ri_yi, yi_pj, 1, 1)
            else
                mul!(ri_pj, ri_yi, yi_pj)
            end
        end
    end

    return J
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian!(J, ::Constant, ::InPlace,
    model::CoupledModel, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydp_cache = similar(J, number_of_inputs(model), number_of_parameters(model))
    )

    # get state, input, and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input jacobian
    get_coupling_parameter_jacobian!(dydp_cache, model, p)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, submodels[i], pi)
        for j = 1:length(submodels)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, submodels[i], pi)
                mul!(ri_pj, ri_yi, yi_pj, 1, 1)
            else
                mul!(ri_pj, ri_yi, yi_pj)
            end
        end
    end

    return J
end

# calculate parameter jacobian for a combination of model
function _get_parameter_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    model::CoupledModel, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(model), number_of_inputs(model)),
    dydp_cache = similar(J, number_of_inputs(model), number_of_parameters(model))
    )

    # get state, input, and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input jacobian
    get_coupling_parameter_jacobian!(dydp_cache, model, dx, x, p, t)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, submodels[i], dxi, xi, yi, pi, t)
        for j = 1:length(submodels)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, submodels[i], dxi, xi, yi, pi, t)
                mul!(ri_pj, ri_yi, yi_pj, 1, 1)
            else
                mul!(ri_pj, ri_yi, yi_pj)
            end
        end
    end

    return J
end

"""
    get_time_gradient(model)
    get_time_gradient(model, p)
    get_time_gradient(model, dx, x, y, p, t)

Calculate the derivative of the residual function for `model` with respect to time
"""
function get_time_gradient(model, args...)
    return _get_time_gradient(inplaceness(model), model, args...)
end

# dispatch to an in-place function
function _get_time_gradient(::InPlace, model)
    Nx = number_of_states(model)
    dT = Vector{Float64}(undef, Nx)
    return get_time_gradient!(dT, model)
end

# dispatch to an in-place function
function _get_time_gradient(::InPlace, model, p)
    TF = eltype(p)
    Nx = number_of_states(model)
    dT = Vector{TF}(undef, Nx)
    return get_time_gradient!(dT, model, p)
end

# dispatch to an in-place function
function _get_time_gradient(::InPlace, model, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(model)
    dT = Vector{TF}(undef, Nx)
    return get_time_gradient!(dT, model, dx, x, y, p, t)
end

# calculate time gradient for a standalone model
function _get_time_gradient(::OutOfPlace, model::Model, args...)
    return gradient(model.tgrad, model.nx, args...)
end

# calculate time derivative for a combination of models
function _get_time_gradient(::OutOfPlace, model::CoupledModel)

    # initialize mass matrix
    r_t = get_time_partial(model)

    # calculate input jacobian
    r_y = get_input_jacobian(model)

    # calculate input jacobian
    y_t = get_coupling_time_gradient(model)

    return r_t + r_y*y_t
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient(::OutOfPlace, model::CoupledModel, p)

    # initialize parameter jacobian
    r_t = get_time_partial(model, p)

    # calculate input jacobian
    r_y = get_input_jacobian(model, p)

    # calculate coupling function parameter jacobian
    y_t = get_coupling_time_gradient(model, p)

    return r_t + r_y*y_t
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient(::OutOfPlace, model::CoupledModel, dx, x, y, p, t)

    # initialize parameter jacobian
    r_t = get_time_partial(model, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(model, dx, x, y, p, t)

    # calculate coupling function parameter jacobian
    y_t = get_coupling_time_gradient(model, dx, x, p, t)

    return r_t + r_y*y_t
end

function get_time_partial(model::CoupledModel) where N
    return vcat(get_time_gradient.(model.models)...)
end

function get_time_partial(model::CoupledModel, p) where N

    ip = parameter_indices(model)

    ps = view.(Ref(ps), ip[i])

    return vcat(get_time_gradient.(model.models, ps)...)
end

function get_time_partial(model::CoupledModel, dx, x, y, p, t) where N

    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    dxs = view.(Ref(dx), ix[i])
    xs = view.(Ref(xs), ix[i])
    ys = view.(Ref(ys), iy[i])
    ps = view.(Ref(ps), ip[i])

    return vcat(get_time_gradient.(model.models, dxs, xs, ys, ps, t)...)
end

"""
    get_time_gradient!(dT, model)
    get_time_gradient!(dT, model, p)
    get_time_gradient!(dT, model, dx, x, y, p, t)

In-place version of [`get_time_gradient`](@ref)
"""
function get_time_gradient!(dT, model, args...; kwargs...)
    return _get_time_gradient!(dT, inplaceness(model), model; kwargs...)
end

# dispatch to an out-of-place function
function _get_time_gradient!(dT, ::OutOfPlace, model, args...; kwargs...)
    return dT .= get_time_gradient(model, args...; kwargs...)
end

# calculate time gradient for a standalone model
function _get_time_gradient!(dT, ::InPlace, model::Model, args...)
    return gradient!(dT, model.tgrad, args...)
end

# calculate parameter jacobian for a coupled model
function _get_time_gradient!(dT, ::InPlace, model::CoupledModel;
    drdy_cache = similar(dT, number_of_states(model), number_of_inputs(model)),
    dydt_cache = similar(dT, number_of_inputs(model))
    )

    # get state and input indices
    ix = state_indices(model)
    iy = input_indices(model)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, model)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, submodels[i])
        for j = 1:length(submodels)
            ri_t = view(dT, ix[i])
            yi_t = view(dydt_cache, iy[i])
            if i == j
                get_parameter_jacobian!(Jij, submodels[i])
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

# calculate parameter jacobian for a coupled model
function _get_time_gradient!(dT, ::InPlace, model::CoupledModel, p;
    drdy_cache = similar(dT, number_of_states(model), number_of_inputs(model)),
    dydt_cache = similar(dT, number_of_inputs(model))
    )

    # get state, input, and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, model, p)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, submodels[i], pi)
        for j = 1:length(submodels)
            ri_t = view(dT, ix[i], ix[j])
            yi_t = view(dydt_cache, iy[i], ix[j])
            if i == j
                get_parameter_jacobian!(ri_t, submodels[i], pi)
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

# calculate parameter jacobian for a coupled model
function _get_time_gradient!(dT, ::InPlace, model::CoupledModel, dx, x, y, p, t;
    drdy_cache = similar(dT, number_of_states(model), number_of_inputs(model)),
    dydt_cache = similar(dT, number_of_inputs(model))
    )

    # get state, input, and parameter indices
    ix = state_indices(model)
    iy = input_indices(model)
    ip = parameter_indices(model)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, model, dx, x, p, t)

    # calculate jacobian
    submodels = model.models
    for i = 1:length(submodels)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, submodels[i], dxi, xi, yi, pi, t)
        for j = 1:length(submodels)
            ri_t = view(dT, ix[i], ix[j])
            yi_t = view(dydt_cache, iy[i], ix[j])
            if i == j
                get_parameter_jacobian!(ri_t, submodels[i], dxi, xi, yi, pi, t)
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

"""
    get_coupling_inputs(model::CoupledModel, dx, x, p, t)

Calculate the coupling inputs for a coupled model.
"""
function get_coupling_inputs(model::CoupledModel, dx, x, p, t)
    get_coupling_inputs(model.coupling, dx, x, p, t)
end

function get_coupling_inputs(coupling::Coupling, dx, x, p, t)
    return _get_coupling_inputs(inplaceness(coupling), coupling, dx, x, p, t)
end

# dispatch to an in-place function
function _get_coupling_inputs(::InPlace, coupling, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(coupling)
    y = zeros(TF, Ny)
    get_coupling_inputs!(y, coupling, dx, x, p, t)
    return y
end

# calculate out-of-place coupling inputs
function _get_coupling_inputs(::OutOfPlace, coupling, dx, x, p, t)
    return coupling.g(dx, x, p, t)
end

"""
    get_coupling_inputs!(y, models::CoupledModel, dx, x, p, t) where N

In-place version of [`get_coupling_inputs`](@ref)
"""
function get_coupling_inputs!(y, model::CoupledModel, dx, x, p, t)
    get_coupling_inputs!(y, model.coupling, dx, x, p, t)
end

function get_coupling_inputs!(y, coupling::Coupling, dx, x, p, t)
    return _get_coupling_inputs!(y, inplaceness(coupling), coupling, dx, x, p, t)
end

# dispatch to an out-of-place function
function _get_coupling_inputs!(y, ::OutOfPlace, coupling, dx, x, p, t)
    return y .= get_coupling_inputs(coupling, dx, x, p, t)
end

# calculate in-place coupling inputs
function _get_coupling_inputs!(y, ::InPlace, coupling, dx, x, p, t) where N
    return coupling.g(y, dx, x, p, t)
end

"""
    get_coupling_rate_jacobian(model)
    get_coupling_rate_jacobian(model, p)
    get_coupling_rate_jacobian(model, dx, x, p, t)

Calculate the jacobian of the coupling function with respect to the state rates
"""
function get_coupling_rate_jacobian(model::CoupledModel, args...)
    return get_coupling_rate_jacobian(model.coupling, args...)
end

function get_coupling_rate_jacobian(coupling::Coupling, args...)
    return _get_coupling_rate_jacobian(inplaceness(coupling), coupling, args...)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::InPlace, coupling)
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{Float64}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, coupling)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::InPlace, coupling, p)
    TF = eltype(p)
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, coupling, p)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::InPlace, coupling, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, coupling, dx, x, p, t)
end

# calculate out-of-place coupling function rate jacobian
function _get_coupling_rate_jacobian(::OutOfPlace, coupling, args...)

    return jacobian(coupling.ratejac, coupling.ny, coupling.nx, args...)
end

"""
    get_coupling_rate_jacobian!(J, model)
    get_coupling_rate_jacobian!(J, model, p)
    get_coupling_rate_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_rate_jacobian`](@ref).
"""
function get_coupling_rate_jacobian!(J, model::CoupledModel, args...)
    return get_coupling_rate_jacobian!(J, model.coupling, args...)
end

function get_coupling_rate_jacobian!(J, coupling::Coupling, args...)
    return _get_coupling_rate_jacobian!(J, inplaceness(coupling), coupling, args...)
end

# dispatch to an out-of-place function
function _get_coupling_rate_jacobian!(J, ::OutOfPlace, coupling, args...)
    return J .= get_coupling_rate_jacobian(coupling, args...)
end

# calculate in-place coupling function rate jacobian
function _get_coupling_rate_jacobian!(J, ::InPlace, coupling, args...)

    return jacobian!(J, coupling.ratejac, args...)
end

"""
    get_coupling_state_jacobian(model)
    get_coupling_state_jacobian(model, p)
    get_coupling_state_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the coupling function for `model` with respect to the
state variables
"""
function get_coupling_state_jacobian(model::CoupledModel, args...)
    return get_coupling_state_jacobian(model.coupling, args...)
end

function get_coupling_state_jacobian(coupling::Coupling, args...)
    return _get_coupling_state_jacobian(inplaceness(coupling), coupling, args...)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::InPlace, coupling)
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{Float64}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, coupling)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::InPlace, coupling, p)
    TF = eltype(p)
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, coupling, p)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::InPlace, coupling, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(coupling)
    Nx = number_of_states(coupling)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, coupling, dx, x, p, t)
end

# calculate out-of-place coupling function state jacobian
function _get_coupling_state_jacobian(::OutOfPlace, coupling, args...)

    return jacobian(coupling.statejac, coupling.ny, coupling.nx, args...)
end

"""
    get_coupling_state_jacobian!(J, model)
    get_coupling_state_jacobian!(J, model, p)
    get_coupling_state_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_state_jacobian`](@ref)
"""
function get_coupling_state_jacobian!(J, model::CoupledModel, args...)
    return get_coupling_state_jacobian!(J, model.coupling, args...)
end

function get_coupling_state_jacobian!(J, coupling, args...)
    return _get_coupling_state_jacobian!(J, inplaceness(coupling), coupling, args...)
end

# dispatch to an out-of-place function
function _get_coupling_state_jacobian!(J, ::OutOfPlace, coupling, args...)
    return J .= get_coupling_state_jacobian(models, args...)
end

# calculate in-place coupling function state jacobian
function _get_coupling_state_jacobian!(J, ::InPlace, coupling, args...)

    return jacobian!(J, coupling.statejac, args...)
end

"""
    get_coupling_parameter_jacobian(model)
    get_coupling_parameter_jacobian(model, p)
    get_coupling_parameter_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the coupling function for `model` with respect to the
parameters
"""
function get_coupling_parameter_jacobian(model::CoupledModel, args...)
    return get_coupling_parameter_jacobian(model.coupling, args...)
end

function get_coupling_parameter_jacobian(coupling::Coupling, args...)
    return _get_coupling_parameter_jacobian(inplaceness(coupling), coupling, args...)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::InPlace, coupling)
    Ny = number_of_inputs(coupling)
    Np = number_of_parameters(coupling)
    J = Matrix{Float64}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, coupling)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::InPlace, coupling, p)
    TF = eltype(p)
    Ny = number_of_inputs(coupling)
    Np = number_of_parameters(coupling)
    J = Matrix{TF}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, coupling)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::InPlace, coupling, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(coupling)
    Np = number_of_parameters(coupling)
    J = Matrix{TF}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, coupling, dx, x, p, t)
end

# calculate out-of-place coupling function parameter jacobian
function _get_coupling_parameter_jacobian(::OutOfPlace, coupling, args...)

    return jacobian(coupling.paramjac, coupling.ny, coupling.np, args...)
end

"""
    get_coupling_parameter_jacobian!(J, model)
    get_coupling_parameter_jacobian!(J, model, p)
    get_coupling_parameter_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_parameter_jacobian`](@ref)
"""
function get_coupling_parameter_jacobian!(J, model::CoupledModel, args...)
    return get_coupling_parameter_jacobian!(J, model.coupling, args...)
end

function get_coupling_parameter_jacobian!(J, coupling::Coupling, args...)
    return _get_coupling_parameter_jacobian!(J, inplaceness(coupling), coupling, args...)
end

# dispatch to an out-of-place function
function _get_coupling_parameter_jacobian!(J, ::OutOfPlace, coupling, args...)
    return J .= get_coupling_parameter_jacobian(coupling, args...)
end

# calculate in-place coupling function parameter jacobian
function _get_coupling_parameter_jacobian!(J, ::InPlace, coupling, args...)

    return jacobian!(J, coupling.paramjac, args...)
end

"""
    get_coupling_time_gradient(model)
    get_coupling_time_gradient(model, p)
    get_coupling_time_gradient(model, dx, x, y, p, t)

Calculate the derivative of the coupling function for `model` with respect to time
"""
function get_coupling_time_gradient(model::CoupledModel, args...)
    return get_coupling_time_gradient(model.coupling, args...)
end

function get_coupling_time_gradient(coupling::Coupling, args...)
    return _get_coupling_time_gradient(inplaceness(coupling), coupling, args...)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::InPlace, coupling)
    Ny = number_of_inputs(coupling)
    dT = Vector{Float64}(undef, Ny)
    return get_coupling_time_gradient!(dT, coupling)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::InPlace, coupling, p)
    TF = eltype(p)
    Ny = number_of_inputs(coupling)
    dT = Vector{TF}(undef, Ny)
    return get_coupling_time_gradient!(dT, coupling, p)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::InPlace, coupling, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(coupling)
    dT = Vector{TF}(undef, Ny)
    return get_coupling_time_gradient!(dT, coupling, dx, x, p, t)
end

# calculate out-of-place coupling function time gradient
function _get_coupling_time_gradient(::OutOfPlace, coupling, args...)

    return gradient(coupling.tgrad, coupling.ny, args...)
end

"""
    get_coupling_time_gradient!(dT, model)
    get_coupling_time_gradient!(dT, model, p)
    get_coupling_time_gradient!(dT, model, dx, x, p, t)

In-place version of [`get_coupling_time_gradient`](@ref)
"""
function get_coupling_time_gradient!(dT, model::CoupledModel, args...)
    return get_coupling_time_gradient!(dT, model.coupling, args...)
end

function get_coupling_time_gradient!(dT, coupling::Coupling, args...)
    return _get_coupling_time_gradient!(dT, inplaceness(coupling), coupling, args...)
end

# dispatch to an out-of-place function
function _get_coupling_time_gradient!(dT, ::OutOfPlace, coupling, args...)
    return dT .= get_coupling_time_gradient(coupling, args...)
end

# calculate in-place coupling function time gradient
function _get_coupling_time_gradient!(dT, ::InPlace, coupling, args...)

    return gradient!(dT, coupling.tgrad, args...)
end

"""
    linearize(model, x, p; kwargs...)

Calculate the jacobian of the residual function for coupled model `model` with respect to 
the state variables and their rates.
"""
function linearize(model::CoupledModel, x, p; dx = FillArrays.Zeros(x), t = 0,
    y = get_coupling_inputs(model, dx, x, p, t))

    return _linearize(inplaceness(model), model, dx, x, y, p, t)
end

function _linearize(::OutOfPlace, model, args...)
    K = Array(get_state_jacobian(model, args...))
    M = Array(get_rate_jacobian(model, args...))
    return K, M
end

function _linearize(::InPlace, model, args...)
    K = get_state_jacobian(model, args...)
    M = get_rate_jacobian(model, args...)
    return K, M
end

"""
    get_eigen(model::TM, K, M; kwargs...)

Return the eigenvalues, left eigenvector matrix, and right eigenvector matrix 
corresponding to the model.

For in-place models, the number of eigenvalues to compute may be specified using
the `nev` keyword argument.
"""
function get_eigen(model, K, M; kwargs...)

    return _get_eigen(inplaceness(model), K, M; kwargs...)
end

function _get_eigen(::OutOfPlace, K, M)
    A = -K
    B = M
    E = eigen(A, B) # eigenvalue decomposition
    λ = eigvals(E) # eigenvalues
    V = eigvecs(E) # right eigenvector matrix
    U = I/(A*V) # left eigenvector matrix
    return λ, U, V
end

function _get_eigen(::InPlace, K, M; nev=min(20, number_of_states(model)))

    # calculate state and rate jacobians
    A = -K
    B = M

    # construct linear map
    nx = size(A, 1)
    TF = promote_type(eltype(A), eltype(B))
    Afact = lu(A)
    f! = (b, x) -> ldiv!(b, Afact, B*x)
    fc! = (b, x) -> mul!(b, B', Afact'\x)
    Abar = LinearMap{TF}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V = partialeigen(partialschur(Abar; nev=min(nx,nev), which=LM())[1])

    # sort eigenvalues by magnitude
    perm = sortperm(λ, by=(λ)->(abs(λ),imag(λ)), rev=true)
    λ = λ[perm[1:nev]]
    V = V[:,perm[1:nev]]

    # eigenvalues are actually 1/λ, no modification necessary for eigenvectors
    λ .= 1 ./ λ

    # also return left eigenvectors
    U = GXBeam.left_eigenvectors(A, -B, λ, V)

    return λ, U, V
end

function correlate_eigenmodes(U, M, V; tracked_modes=1:size(U,1), rtol=1e-3, atol=1e-3)

    # get number of (not necessarily unique) eigenvalues
    nev = size(U,1)

    # construct the correlation matrix
    C = U*M*V

    # get the row permutation that puts maximum values on the diagonals
    rval = Vector{real(eltype(C))}(undef, nev) # magnitude of values in chosen row
    cval = Vector{real(eltype(C))}(undef, nev) # magnitude of values in chosen column
    perm = zeros(Int, nev) # permutation
    corruption = Vector{real(eltype(C))}(undef, nev) # corruption index for each correlation

    # tracked modes
    for i in tracked_modes

        # get correlation matrix entry magnitudes for this row
        for j = 1:nev
            rval[j] = abs(C[i,j])
        end

        # rank modes correlations from best to worst for this row
        ranked_modes = sortperm(rval, rev=true)

        # choose the best unassigned mode correlation for this row
        j1 = ranked_modes[findfirst((x) -> !(x in perm), ranked_modes)]

        # get correlation matrix entry magnitudes for the chosen column
        for k = 1:nev
            cval[k] = abs(C[k,j1])
        end

        # find approximate multiplicity of left and/or right eigenvector
        multiplicity = max(
            count(x -> isapprox(x, rval[j1]; rtol, atol), rval),
            count(x -> isapprox(x, cval[i]; rtol, atol), cval)
        )

        if multiplicity < nev
            # determine n+1 best unassigned fit, where n is the approximate multiplicity
            if j1 in ranked_modes[1:multiplicity]
                j2 = ranked_modes[multiplicity+1]
            else
                j2 = ranked_modes[1]
            end

            # assign best eigenmode fit, create corruption index
            perm[i] = j1
            corruption[i] = rval[j2]/rval[j1]
        else
            # assign best eigenmode fit, create corruption index
            perm[i] = j1
            corruption[i] = 0.0
        end
    end

    # remaining modes
    for i = 1:nev

        if !(i in tracked_modes)

            # get correlation matrix entry magnitudes for this row
            for j = 1:nev
                rval[j] = abs(C[i,j])
            end

            # rank modes correlations from best to worst for this row
            ranked_modes = sortperm(rval, rev=true)

            # choose the best unassigned mode correlation for this row
            j1 = ranked_modes[findfirst((x) -> !(x in perm), ranked_modes)]

            # get correlation matrix entry magnitudes for the chosen column
            for k = 1:nev
                cval[k] = abs(C[k,j1])
            end

            # find approximate multiplicity of left and/or right eigenvector
            multiplicity = max(
                count(x -> isapprox(x, rval[j1]; rtol, atol), rval),
                count(x -> isapprox(x, cval[i]; rtol, atol), cval)
            )

            # determine n+1 best unassigned fit, where n is the approximate multiplicity
            if j1 in ranked_modes[1:multiplicity]
                j2 = ranked_modes[multiplicity+1]
            else
                j2 = ranked_modes[1]
            end

            # assign best eigenmode fit, create corruption index
            perm[i] = j1
            corruption[i] = rval[j2]/rval[j1]
        end

    end

    return perm, corruption
end


"""
    get_ode(model)
    get_ode(model, p)

Construct an ODEFunction corresponding to the specified model or models which
may be solved using DifferentialEquations.
"""
function get_ode(model::CoupledModel, p=nothing)

    iip = isinplace(model)
    compile = true

    # construct coupling function
    fy = ode_coupling_function(model)

    # construct rate function
    f = ode_rate_function(model, fy)

    # construct mass matrix
    mass_matrix = ode_mass_matrix(model, fy)

    # construct time gradient function
    tgrad = ode_time_gradient(model, fy)
 
    # construct jacobian function
    jac = ode_state_jacobian(model, fy)

    # construct parameter jacobian function
    paramjac = ode_parameter_jacobian(model, fy)

    # return ODE function
    return ODEFunction{iip, compile}(f; mass_matrix, tgrad, jac, paramjac)
end

function ode_coupling_function(model::CoupledModel)

    # check whether cached variables should be used
    if isinplace(model)
        # model is inplace so use cached variables if possible

        # problem dimensions
        Nx = number_of_states(model)
        Ny = number_of_inputs(model)
        Np = number_of_parameters(model)

        # cached variables
        ucache = fill(NaN, Nx)
        ycache = fill(NaN, Ny)
        pcache = fill(NaN, Np)
        tcache = fill(NaN)

        # state rate vector
        dx = FillArrays.Zeros(x)

        # coupling function
        fy = (x, p, t) -> begin
            # check if we can use the cache variables (no custom types)
            if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
                # update the cache variables (if necessary)
                if (x != ucache) && (p != pcache) && (t != tcache[])
                    ucache .= x
                    ycache .= get_coupling_inputs!(ycache, model, dx, x, p, t)
                    pcache .= p
                    tcache .= t
                end
                # set inputs to the cached inputs
                y = ycache
            else
                # calculate inputs (out-of-place to accomodate custom type)
                y = get_coupling_inputs(model, dx, x, p, t)
            end
            # return cached or computed inputs
            return y
        end
    else
        # model is out of place so calculate coupling inputs directly
        fy = (x, p, t) -> get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
    end

    # return coupling function
    return fy
end

function ode_rate_function(model::CoupledModel, fy)

    if isinplace(model)
        f = (dx, x, p, t) -> begin
            dx = FillArrays.Zeros(x)
            y = fy(x, p, t)
            get_residual!(resid, model, dx, x, y, p, t)
            return resid .*= -1
        end
    else
        f = (x, p, t) -> begin
            dx = FillArrays.Zeros(x)
            y = fy(x, p, t)
            resid = get_residual(model, dx, x, y, p, t)
            return -resid
        end
    end

    return f
end

function ode_mass_matrix(model::CoupledModel, fy, p=nothing)

    submodels = model.models
    coupling = model.coupling

    if all(isidentity.(getproperty.(submodels, :ratejac))) && iszero(coupling.ratejac)
        # mass matrix is the identity matrix
        mass_matrix = I
    elseif all(isinvariant.(getproperty.(submodels, :ratejac))) && isinvariant(coupling.ratejac)
        # mass matrix is independent of all inputs
        mass_matrix = get_rate_jacobian(model)
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :ratejac))) && 
        isconstant(coupling.ratejac)
        # mass matrix only depends on the parameter vector
        mass_matrix = get_rate_jacobian(model, p)
    else
        # mass matrix defined as a function of state variables, parameters, and time 

        # problem dimensions
        Nx = number_of_states(model)
        Ny = number_of_inputs(model)

        # initialize mass matrix
        M = zeros(Nx, Nx)

        # construct update function
        if isinplace(model)
            # model is inplace so use cached variables if possible

            # cached variables
            drdy_cache = fill(NaN, Nx, Ny),
            dyddx_cache = fill(NaN, Ny, Nx)

            # state rate vector
            dx = FillArrays.Zeros(x)

            # update function
            update_func = (M, x, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
                    # use cache variables
                    get_rate_jacobian!(M, model, dx, x, fy(x, p, t), p, t;
                        drdy_cache = drdy_cache, dyddx_cache = dyddx_cache)
                else
                    # don't use cache variables (to accomodate custom types)
                    get_rate_jacobian!(M, model, dx, x, fy(x, p, t), p, t)
                end
            end
        else
            # model is out-of-place so don't use cache variables
            update_func = (M, x, p, t) -> M .= get_rate_jacobian(model, dx, x, fy(x, p, t), p, t)
        end

        mass_matrix = DiffEqArrayOperator(M; update_func)
    end

    return mass_matrix 
end

function ode_time_gradient(model::CoupledModel, fy, p=nothing)

    submodels = model.models
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :tgrad))) && iszero(coupling.tgrad)
        # time gradient is a zero vector
        if isinplace(model)
            tgrad = (dT, x, p, t) -> dT .= 0
        else
            tgrad = (x, p, t) -> zero(x)
        end
    elseif all(isinvariant.(getproperty.(submodels, :tgrad))) && isinvariant(coupling.tgrad)
        # time gradient is independent of all inputs
        dT0 = get_time_gradient(model)
        if isinplace(model)
            tgrad = (dT, x, p, t) -> dT .= dT0
        else
            tgrad = (x, p, t) -> dT0
        end
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :tgrad))) && 
        isconstant(coupling.tgrad)
        # time gradient only depends on the parameter vector
        dT0 = get_rate_jacobian(model, p)
        if isinplace(model)
            tgrad = (dT, x, p, t) -> dT .= dT0
        else
            tgrad = (x, p, t) -> dT0
        end
    else
        # time gradient defined as a function of state variables, parameters, and time 

        if isinplace(model)
            # model is inplace so use cached variables if possible

            # problem dimensions
            Nx = number_of_states(model)
            Ny = number_of_inputs(model)

            # cached variables
            drdy_cache = fill(NaN, Nx, Ny),
            dyddx_cache = fill(NaN, Ny, Nx)

            # state rate vector
            dx = FillArrays.Zeros(x)

            # update function
            tgrad = (dT, x, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
                    # use cache variables
                    get_time_gradient!(dT, model, dx, x, fy(x, p, t), p, t;
                        drdy_cache = drdy_cache, dyddx_cache = dyddx_cache)
                else
                    # don't use cache variables (to accomodate custom types)
                    get_time_gradient!(dT, model, dx, x, fy(x, p, t), p, t)
                end
            end
        else
            # model is out-of-place so don't use cache variables
            tgrad = (x, p, t) -> get_time_gradient(model, dx, x, fy(x, p, t), p, t)
        end
    end

    return tgrad
end

function ode_state_jacobian(model::CoupledModel, fy, p=nothing)

    submodels = model.models
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :statejac))) && iszero(coupling.statejac)
        # state jacobian is a zero vector
        if isinplace(model)
            jac = (J, x, p, t) -> J .= 0
        else
            jac = (x, p, t) -> zero(x) .* zero(x)'
        end
    elseif all(isinvariant.(getproperty.(submodels, :statejac))) && isinvariant(coupling.statejac)
        # state jacobian is independent of all inputs
        J0 = get_state_jacobian(model)
        if isinplace(model)
            jac = (J, x, p, t) -> J .= J0
        else
            jac = (x, p, t) -> J0
        end
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :statejac))) && 
        isconstant(coupling.statejac)
        # state jacobian only depends on the parameter vector
        J0 = get_state_jacobian(model, p)
        if isinplace(model)
            jac = (J, x, p, t) -> J .= J0
        else
            jac = (x, p, t) -> J0
        end
    else
        # state jacobian defined as a function of state variables, parameters, and time 

        # construct state jacobian function
        if isinplace(model)
            # model is inplace so use cached variables if possible

            # problem dimensions
            Nx = number_of_states(model)
            Ny = number_of_inputs(model)

            # cached variables
            drdy_cache = fill(NaN, Nx, Ny),
            dyddx_cache = fill(NaN, Ny, Nx)

            # state rate vector
            dx = FillArrays.Zeros(x)

            # update function
            jac = (J, x, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
                    # use cache variables
                    get_state_jacobian!(J, model, dx, x, fy(x, p, t), p, t;
                        drdy_cache = drdy_cache, dyddx_cache = dyddx_cache)
                else
                    # don't use cache variables (to accomodate custom types)
                    get_state_jacobian!(J, model, dx, x, fy(x, p, t), p, t)
                end
            end
        else
            # model is out-of-place so don't use cache variables
            jac = (x, p, t) -> get_state_jacobian(model, dx, x, fy(x, p, t), p, t)
        end
    end

    return jac
end

function ode_parameter_jacobian(model::CoupledModel, fy, p=nothing)

    submodels = model.models
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :paramjac))) && isidentity(coupling.paramjac)
        # parameter jacobian is a zero matrix
        if isinplace(model)
            jac = (pJ, x, p, t) -> pJ .= 0
        else
            jac = (x, p, t) -> zero(x) .* zero(p)'
        end
    elseif all(isinvariant.(getproperty.(submodels, :paramjac))) && isinvariant(coupling.paramjac)
        # parameter jacobian is independent of all inputs
        pJ0 = get_parameter_jacobian(model)
        if isinplace(model)
            paramjac = (pJ, x, p, t) -> pJ .= pJ0
        else
            paramjac = (x, p, t) -> pJ0
        end
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :paramjac))) && 
        isconstant(coupling.paramjac)
        # time gradient only depends on the parameter vector
        pJ0 = get_parameter_jacobian(model, p)
        if isinplace(model)
            paramjac = (pJ, x, p, t) -> pJ .= pJ0
        else
            paramjac = (x, p, t) -> pJ0
        end
    else
        # construct time gradient function
        if isinplace(model)
            # model is inplace so use cached variables if possible

            # problem dimensions
            Nx = number_of_states(model)
            Ny = number_of_inputs(model)

            # cached variables
            drdy_cache = fill(NaN, Nx, Ny),
            dyddx_cache = fill(NaN, Ny, Nx)

            # state rate vector
            dx = FillArrays.Zeros(x)

            # update function
            pjac = (dT, x, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
                    # use cache variables
                    get_parameter_jacobian!(pJ, model, dx, x, fy(x, p, t), p, t;
                        drdy_cache = drdy_cache, dyddx_cache = dyddx_cache)
                else
                    # don't use cache variables (to accomodate custom types)
                    get_parameter_jacobian!(pJ, model, dx, x, fy(x, p, t), p, t)
                end
            end
        else
            # model is out-of-place so don't use cache variables
            pjac = (x, p, t) -> get_parameter_jacobian(model, dx, x, fy(x, p, t), p, t)
        end
    end

    return pjac
end