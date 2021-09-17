"""
    couple_models(models...)

Couples multiple models together to form a coupled model.  The coupled model is
represented as a tuple of models, though this may change in the future.
"""
couple_models(models...) = models

"""
    number_of_states(model)

Return the total number of state variables corresponding to the model or models.
"""
number_of_states

# by default, dispatch on type
number_of_states(model::TM) where TM <: AbstractModel = number_of_states(TM)

# models with no states have... no states
number_of_states(::Type{TM}) where TM <: NoStateModel = 0

# coupled models have concatenated states
function number_of_states(models::NTuple{N,AbstractModel}) where N
    return sum(number_of_states.(models))
end

# coupled models have concatenated states
function number_of_states(::Type{TM}) where TM <: NTuple{N,AbstractModel} where N
    return sum(number_of_states.(TM.parameters))
end

"""
    number_of_inputs(model)

Return the total number of inputs corresponding to the model or models.
"""
number_of_inputs

# by default, dispatch on type
number_of_inputs(model::TM) where TM <: AbstractModel = number_of_inputs(TM)

# models with no state variables also have no inputs
number_of_inputs(::Type{TM}) where TM <: NoStateModel = 0

# coupled models have concatenated inputs
function number_of_inputs(models::NTuple{N,AbstractModel}) where N
    sum(number_of_inputs.(models))
end

# coupled models have concatenated inputs
function number_of_inputs(::Type{TM}) where TM <: NTuple{N,AbstractModel} where N
    sum(number_of_inputs.(TM.parameters))
end

"""
    number_of_parameters(model)

Return the total number of parameters corresponding to the model or models.
"""
number_of_parameters

# by default, dispatch on type
number_of_parameters(model::TM) where TM <: AbstractModel = number_of_parameters(TM)

# coupled models have concatenated parameters + additional parameters
function number_of_parameters(models::NTuple{N,AbstractModel}) where N
    sum(number_of_parameters.(models)) + number_of_additional_parameters(models...)
end

# coupled models have concatenated parameters + additional parameters
function number_of_parameters(::Type{TM}) where TM <: NTuple{N,AbstractModel} where N
    return sum(number_of_parameters.(TM.parameters)) + number_of_additional_parameters(TM.parameters...)
end

"""
    number_of_additional_parameters(models...)

Return the total number of additional parameters corresponding to the model coupling.
"""
number_of_additional_parameters

function number_of_additional_parameters(models::Vararg{AbstractModel,N}) where N
    return number_of_additional_parameters(typeof.(models)...)
end

function number_of_additional_parameters(models::Type{TM}) where TM <: NTuple{N,AbstractModel} where N
    return number_of_parameters(TM.parameters...)
end

function number_of_additional_parameters(models::NTuple{N,AbstractModel}) where N
    return number_of_parameters(models...)
end

"""
    state_indices(models)

Return the indices corresponding to the state variables for each model in `models`
"""
function state_indices(models::NTuple{N,AbstractModel}) where N
    # input indices
    Nx = number_of_states.(models)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    return UnitRange.(ix1, ix2)
end

"""
    input_indices(models)

Return the indices corresponding to the input variables for each model in `models`
"""
function input_indices(models::NTuple{N,AbstractModel}) where N
    # input indices
    Ny = number_of_inputs.(models)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    return UnitRange.(iy1, iy2)
end

"""
    parameter_indices(models)

Return the indices corresponding to the parameters for each model in `models`
"""
function parameter_indices(models::NTuple{N,AbstractModel}) where N
    # parameter indices
    Np = number_of_parameters.(models)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    return UnitRange.(ip1, ip2)
end

"""
    additional_parameter_indices(models)

Return the indices of the additional parameters for the coupling function
corresponding to the models in `models`
"""
function additional_parameter_indices(models)
    Npi = number_of_parameters.(models)
    Np = number_of_parameters(models)
    ip1 = sum(Npi) + 1
    ip2 = Np
    return UnitRange(ip1, ip2)
end

"""
    get_states(model; kwargs...)

Return the state variable vector corresponding to `model` using the state
variable values in `kwargs`.
"""
function get_states(model::AbstractModel; kwargs...)
    x = zeros(number_of_states(model))
    set_states!(x, model; kwargs...)
    return x
end

"""
    get_inputs(model; kwargs...)

Return the input vector corresponding to `model` using the input values in `kwargs`.
"""
function get_inputs(model::AbstractModel; kwargs...)
    y = zeros(number_of_inputs(model))
    set_inputs!(y, model; kwargs...)
    return y
end

"""
    get_parameters(model; kwargs...)

Return the parameter vector corresponding to `model` using the parameter values
in `kwargs`.
"""
function get_parameters(model::AbstractModel; kwargs...)
    p = zeros(number_of_parameters(model))
    set_parameters!(p, model; kwargs...)
    return p
end

"""
    get_additional_parameters(model; kwargs...)

Return the elements of the parameter vector corresponding to the additional
parameters of coupled model `model` using the parameters specified in `kwargs`
"""
function get_additional_parameters(models::NTuple{N,AbstractModel};
    kwargs...) where N

    padd = zeros(number_of_additional_parameters(models...))

    set_additional_parameters!(padd, models...; kwargs...)

    return padd
end

"""
    set_states!(x, model; kwargs...)

In-place version of [`get_states`](@ref)
"""
set_states!

set_states!(x, model::NoStateModel) = x

"""
    set_states!(x, model, i; kwargs...)

Set the elements of the state variable vector `x` corresponding to the `i`th
model of coupled model `model` to the values specified in `kwargs`
"""
function set_states!(x, model::NTuple{N,AbstractModel}, i; kwargs...) where N

    return set_states!(view(x, state_indices(model)[i]), model[i]; kwargs...)
end

"""
    set_inputs!(y, model; kwargs...)

In-place version of [`get_coupling_inputs`](@ref)
"""
set_inputs!

set_inputs!(y, model::NoStateModel) = y

"""
    set_inputs!(y, model, i; kwargs...)

Set the elements of the input vector `y` corresponding to the `i`th
model of coupled model `model` to the values specified in `kwargs`
"""
function set_inputs!(y, model::NTuple{N,AbstractModel}, i; kwargs...) where N

    return set_inputs!(view(y, input_indices(model)[i]), model[i]; kwargs...)
end

"""
    set_parameters!(p, model; kwargs...)

In-place version of [`get_parameters`](@ref)
"""
set_parameters!

"""
    set_parameters!(p, model, i; kwargs...)

Set the elements of the parameter vector `p` corresponding to the `i`th
model of coupled model `model` to the values specified in `kwargs`
"""
function set_parameters!(p, model::NTuple{N,AbstractModel}, i; kwargs...) where N

    return set_parameters!(view(p, parameter_indices(model)[i]), model[i]; kwargs...)
end

"""
    set_additional_parameters!(p, model; kwargs...)

Set the elements of the parameter vector `p` corresponding to the additional
parameters of coupled model `model` to the values specified in `kwargs`
"""
function set_additional_parameters!(p, models::NTuple{N,AbstractModel};
    kwargs...) where N

    return set_additional_parameters!(view(p, additional_parameter_indices(models)),
        models...; kwargs...)
end

"""
    separate_states(model, x)

Separate the state vector entries in `x` corresponding to each model in `model`
"""
separate_states

separate_states(model::NoStateModel, x) = ()

function separate_states(models::NTuple{N,AbstractModel}, x) where N

    xs = view.(Ref(x), state_indices(models))

    return ntuple(i->separate_states(models[i], xs[i]), N)
end

"""
    separate_inputs(model, y)

Separate the input vector entries in `y` corresponding to each model in `model`
"""
separate_inputs

separate_inputs(model::NoStateModel, y) = ()

function separate_inputs(models::NTuple{N,AbstractModel}, y) where N

    ys = view.(Ref(y), input_indices(models))

    return ntuple(i->separate_inputs(models[i], ys[i]), N)
end

"""
    separate_parameters(model, p)

Separate the parameter vector entries in `p` corresponding to each model in `model`
and the additional parameters corresponding to the model coupling, if applicable.
"""
separate_parameters

function separate_parameters(models::NTuple{N,AbstractModel}, p) where N

    pmodels = view.(Ref(p), parameter_indices(models))
    pmodels_sep = ntuple(i->separate_parameters(models[i], pmodels[i]), N)

    padd = view(p, additional_parameter_indices(models))
    padd_sep = separate_additional_parameters(models..., padd)

    return (pmodels_sep..., padd_sep)
end

"""
    separate_additional_parameters(model, p)

Separate the additional parameter vector entries in `p` corresponding to the
model coupling, if applicable.
"""
separate_additional_parameters

"""
    get_residual(models, dx, x, y, p, t)

Calculate the residual for the specified model or models.
"""
function get_residual(models::T, dx, x, y, p, t) where T
    return _get_residual(inplaceness(T), models, dx, x, y, p, t)
end

# models with no states have no residual
function get_residual(models::T, args...) where T <: NoStateModel
    return SVector{0,Float64}()
end

# dispatch to an in-place function
function _get_residual(::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    resid = similar(dx, TF)
    get_residual!(resid, models, dx, x, y, p, t)
    return resid
end

# calculate residual for a combination of models
function _get_residual(::OutOfPlace, models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    return get_model_residual(models, dx, x, y, p, t)
end

@generated function get_model_residual(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize residual names
    resid = [Symbol("resid", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(resid[i]) = get_residual(models[$i], dx[$(ix[i])], x[$(ix[i])],
                y[$(iy[i])], p[$(ip[i])], t)
        end
    end

    expr = quote
        $expr
        resid = vcat($(resid...))
    end

    return expr
end

"""
    get_residual!(resid, models, dx, x, y, p, t)

In-place version of [`get_residual`](@ref)
"""
function get_residual!(resid, models::T, dx, x, y, p, t) where T
    return _get_residual!(resid, inplaceness(T), models, dx, x, y, p, t)
end

# dispatch to an out-of-place function
function _get_residual!(resid, ::OutOfPlace, models, dx, x, y, p, t)
    return resid .= get_residual(models, dx, x, y, p, t)
end

# calculate residual for a combination of models
function _get_residual!(resid, ::InPlace, models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    for i = 1:length(models)
        vresid = view(resid, ix[i])
        vdx = view(dx, ix[i])
        vx = view(x, ix[i])
        vy = view(y, iy[i])
        vp = view(p, ip[i])
        get_residual!(vresid, models[i], vdx, vx, vy, vp, t)
    end

    return resid
end

# --- Performance Overloads --- #

"""
    get_rate_jacobian(model)
    get_rate_jacobian(model, p)
    get_rate_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the residual function for `model` with respect to the
state rates
"""
function get_rate_jacobian(models::TM, args...) where TM
    return _get_rate_jacobian(rate_jacobian_type(TM), inplaceness(TM), models, args...)
end

# dispatch to an in-place function
function _get_rate_jacobian(::Any, ::InPlace, models)
    Nx = number_of_states(models)
    J = Matrix{Float64}(undef, Nx, Nx)
    return get_rate_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_rate_jacobian(::Any, ::InPlace, models, p)
    TF = promote_type(eltype(p))
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_rate_jacobian!(J, models, p)
end

# dispatch to an in-place function
function _get_rate_jacobian(::Any, ::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_rate_jacobian!(J, models, dx, x, y, p, t)
end

# return an empty matrix
function _get_rate_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return zeros(SMatrix{Nx, Nx, Float64})
end

# return a zero matrix
function _get_rate_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return zeros(SMatrix{Nx, Nx, Float64})
end

# return the identity matrix
function _get_rate_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return SMatrix{Nx, Nx, Float64}(I)
end

# dispatch to `get_rate_jacobian` without arguments
function _get_rate_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_rate_jacobian(models)
end

# dispatch to `get_rate_jacobian` without arguments
function _get_rate_jacobian(::Invariant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_rate_jacobian(models)
end

# dispatch to `get_rate_jacobian` with only parameters as arguments
function _get_rate_jacobian(::Constant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_rate_jacobian(models, p)
end

# use automatic differentiation since a custom jacobian definition is absent
function _get_rate_jacobian(::Invariant, ::OutOfPlace, models::TM) where TM

    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)
    Np = number_of_parameters(TM)

    dx = zeros(SVector{Nx,Float64})
    x = zeros(SVector{Nx,Float64})
    y = zeros(SVector{Ny,Float64})
    p = zeros(SVector{Np,Float64})
    t = 0

    f = dx -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, dx)
end

# use automatic differentiation since a custom jacobian definition is absent
function _get_rate_jacobian(::Constant, ::OutOfPlace, models::TM, p) where TM

    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)

    dx = zeros(SVector{Nx,Float64})
    x = zeros(SVector{Nx,Float64})
    y = zeros(SVector{Ny,Float64})
    t = 0

    f = dx -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, dx)
end

# use automatic differentiation since a custom jacobian definition is absent
function _get_rate_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace, models, dx, x, y, p, t)

    f = dx -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, dx)
end

# calculate the rate jacobian for a combination of models
function _get_rate_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(models)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(models)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(models)

    # calculate jacobian (chain rule)
    return r_dx + r_y*y_dx
end

# calculate the rate jacobian for a combination of models
function _get_rate_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, p) where N

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(models, p)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(models, p)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(models, p)

    # calculate jacobian (chain rule)
    return r_dx + r_y*y_dx
end

# calculate the rate jacobian for a combination of models
function _get_rate_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, dx, x, y, p, t) where N

    # partial of residual expression wrt state rates
    r_dx = get_rate_partial(models, dx, x, y, p, t)

    # partial of residual expression wrt to inputs
    r_y = get_input_jacobian(models, dx, x, y, p, t)

    # partial of input expression wrt state rates
    y_dx = get_coupling_rate_jacobian(models, dx, x, p, t)

    return r_dx + r_y*y_dx
end

function get_rate_partial(models::NTuple{N,AbstractModel}) where N
    get_invariant_rate_partial(models)
end

function get_rate_partial(models::NTuple{N,AbstractModel}, p) where N
    get_constant_rate_partial(models, p)
end

function get_rate_partial(models::NTuple{N,AbstractModel}, dx, x, y, p, t) where N
    get_varying_rate_partial(models, dx, x, y, p, t)
end

@generated function get_invariant_rate_partial(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nx = number_of_states.(models.parameters)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_rate_jacobian(models[1]), $(Jij[1, 2:end]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[i, 1:i-1]...), get_rate_jacobian(models[$i]), $(Jij[i, i+1:end]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_constant_rate_partial(models::NTuple{N,AbstractModel},
    p) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_rate_jacobian(models[1], p[$(ip[1])]), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_rate_jacobian(models[$i],
                p[$(ip[i])]), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_varying_rate_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_rate_jacobian(models[1], dx[$(ix[1])], x[$(ix[1])],
            y[$(iy[1])], p[$(ip[1])], t), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_rate_jacobian(models[$i],
                dx[$(ix[i])], x[$(ix[i])], y[$(iy[i])], p[$(ip[i])], t),
                $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

"""
    get_rate_jacobian!(J, model)
    get_rate_jacobian!(J, model, p)
    get_rate_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_rate_jacobian`](@ref).
"""
function get_rate_jacobian!(J, model::TM; kwargs...) where TM
    return _get_rate_jacobian!(J, rate_jacobian_type(TM), inplaceness(TM),
        model; kwargs...)
end

function get_rate_jacobian!(J, model::TM, p; kwargs...) where TM
    return _get_rate_jacobian!(J, rate_jacobian_type(TM), inplaceness(TM),
        model, p; kwargs...)
end

function get_rate_jacobian!(J, model::TM, dx, x, y, p, t; kwargs...) where TM
    return _get_rate_jacobian!(J, rate_jacobian_type(TM), inplaceness(TM),
        model, dx, x, y, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_rate_jacobian!(J, ::Any, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_rate_jacobian(model, args...; kwargs...)
end

# return an empty matrix
function _get_rate_jacobian!(J, ::Empty, ::InPlace, model, args...; kwargs...)
    return J
end

# return a zero matrix
function _get_rate_jacobian!(J, ::Zeros, ::InPlace, model, args...; kwargs...)
    return J .= 0
end

# return the identity matrix
function _get_rate_jacobian!(J, ::Identity, ::InPlace, model, args...; kwargs...)
    J .= 0
    for i = 1:number_of_states(model)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_rate_jacobian!` without arguments
function _get_rate_jacobian!(J, ::Invariant, ::InPlace, model, dx, x, y, p, t; kwargs...)

    return get_rate_jacobian!(J, model; kwargs...)
end

# dispatch to `get_rate_jacobian!` without arguments
function _get_rate_jacobian!(J, ::Invariant, ::InPlace, model, p; kwargs...)

    return get_rate_jacobian!(J, model; kwargs...)
end

# dispatch to `get_rate_jacobian!` with only parameters as arguments
function _get_rate_jacobian!(J, ::Constant, ::InPlace, model, dx, x, y, p, t; kwargs...)

    return get_rate_jacobian!(J, model, p; kwargs...)
end

# use automatic differentiation since a custom definition is absent
function _get_rate_jacobian!(J, ::Invariant, ::InPlace, model)

    dx = FillArrays.Zeros(number_of_states(model))
    x = FillArrays.Zeros(number_of_states(model))
    y = FillArrays.Zeros(number_of_inputs(model))
    p = FillArrays.Zeros(number_of_parameters(model))
    t = 0

    f = dx -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, dx))
end

# use automatic differentiation since a custom definition is absent
function _get_rate_jacobian!(J, ::Constant, ::InPlace, model, p)

    dx = FillArrays.Zeros(number_of_states(model))
    x = FillArrays.Zeros(number_of_states(model))
    y = FillArrays.Zeros(number_of_inputs(model))
    t = 0

    f = dx -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, dx))
end

# use automatic differentiation since a custom definition is absent
function _get_rate_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace, model, dx, x, y, p, t)

    f = dx -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, dx))
end

# calculate rate jacobian for a combination of models
function _get_rate_jacobian!(J, ::Invariant, ::InPlace,
    models::NTuple{N, AbstractModel};
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dyddx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state and input indices
    ix = state_indices(models)
    iy = input_indices(models)

    # calculate input mass matrix
    get_coupling_rate_jacobian!(dyddx_cache, models)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, models[i])
        for j = 1:length(models)
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, models)
                mul!(ri_dxj, ri_yi, yi_dxj, 1, 1)
            else
                mul!(ri_dxj, ri_yi, yi_dxj)
            end
        end
    end

    return J
end

# calculate rate jacobian for a combination of models
function _get_rate_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N, AbstractModel}, p;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dyddx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input rate jacobian
    get_coupling_rate_jacobian!(dyddx_cache, models, p)

    # calculate rate jacobian
    for i = 1:N
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        get_input_jacobian!(ri_yi, models[i], p)
        for j = 1:N
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, models[i], p)
                mul!(ri_dxj, ri_yi, yi_dxj, 1, 1)
            else
                mul!(ri_dxj, ri_yi, yi_dxj)
            end
        end
    end

    return J
end

# calculate rate jacobian for a combination of models
function _get_rate_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N, AbstractModel}, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dyddx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input rate jacobian
    get_coupling_rate_jacobian!(dyddx_cache, models, dx, x, p, t)

    # calculate rate jacobian
    for i = 1:N
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], dxi, xi, yi, pi, t)
        for j = 1:N
            ri_dxj = view(J, ix[i], ix[j])
            yi_dxj = view(dyddx_cache, iy[i], ix[j])
            if i == j
                get_rate_jacobian!(ri_dxj, models[i], dxi, xi, yi, pi, t)
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
function get_state_jacobian(model::TM, args...) where TM
    return _get_state_jacobian(state_jacobian_type(TM), inplaceness(TM), model, args...)
end

# dispatch to an in-place function
function _get_state_jacobian(::Any, ::InPlace, models)
    Nx = number_of_states(models)
    J = Matrix{Float64}(undef, Nx, Nx)
    return get_state_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_state_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_state_jacobian!(J, models, p)
end

# dispatch to an in-place function
function _get_state_jacobian(::Any, ::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Nx, Nx)
    return get_state_jacobian!(J, models, dx, x, y, p, t)
end

# return an empty matrix
function _get_state_jacobian(::Empty, ::OutOfPlace, models, args...)
    return zeros(SMatrix{0, 0, Float64})
end

# return a zero matrix
function _get_state_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return zeros(SMatrix{Nx, Nx, Float64})
end

# return the identity matrix
function _get_state_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return SMatrix{Nx, Nx, Float64}(I)
end

# dispatch to the `get_state_jacobian` without arguments
function _get_state_jacobian(::Invariant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_state_jacobian(models)
end

# dispatch to the `get_state_jacobian` without arguments
function _get_state_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_state_jacobian(models)
end

# dispatch to `get_state_jacobian` with only parameters as arguments
function _get_state_jacobian(::Constant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_state_jacobian(models, p)
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian(::Invariant, ::OutOfPlace, model::TM) where TM

    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)
    Np = number_of_parameters(TM)

    dx = zeros(SVector{Nx,Float64})
    x = zeros(SVector{Nx,Float64})
    y = zeros(SVector{Ny,Float64})
    p = zeros(SVector{Np,Float64})
    t = 0

    f = x -> get_residual(model, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, x)
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian(::Constant, ::OutOfPlace, models, p)

    f = x -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, x)
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace, models, dx, x, y, p, t)

    f = x -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, x)
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    r_x = get_state_partial(models)

    # calculate input jacobian
    r_y = get_input_jacobian(models)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(models)

    return r_x + r_y*y_x
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, p) where N

    # initialize mass matrix
    r_x = get_state_partial(models, p)

    # calculate input jacobian
    r_y = get_input_jacobian(models, p)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(models, p)

    return r_x + r_y*y_x
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, dx, x, y, p, t) where N

    # initialize mass matrix
    r_x = get_state_partial(models, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(models, dx, x, y, p, t)

    # calculate input jacobian
    y_x = get_coupling_state_jacobian(models, dx, x, p, t)

    return r_x + r_y*y_x
end

function get_state_partial(models::NTuple{N,AbstractModel}) where N
    get_invariant_state_partial(models)
end

function get_state_partial(models::NTuple{N,AbstractModel}, p) where N
    get_constant_state_partial(models, p)
end

function get_state_partial(models::NTuple{N,AbstractModel}, dx, x, y, p, t) where N
    get_varying_state_partial(models, dx, x, y, p, t)
end

@generated function get_invariant_state_partial(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nx = number_of_states.(models.parameters)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_state_jacobian(models[1]), $(Jij[1, 2:end]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[i, 1:i-1]...), get_state_jacobian(models[$i]), $(Jij[i, i+1:end]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_constant_state_partial(models::NTuple{N,AbstractModel},
    p) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_state_jacobian(models[1], p[$(ip[1])]), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_state_jacobian(models[$i],
                p[$(ip[i])]), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_varying_state_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Nx[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_state_jacobian(models[1], dx[$(ix[1])], x[$(ix[1])],
            y[$(iy[1])], p[$(ip[1])], t), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_state_jacobian(models[$i],
                dx[$(ix[i])], x[$(ix[i])], y[$(iy[i])], p[$(ip[i])], t), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

"""
    get_state_jacobian!(J, model)
    get_state_jacobian!(J, model, p)
    get_state_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_state_jacobian`](@ref)
"""
function get_state_jacobian!(J, model::TM; kwargs...) where TM
    return _get_state_jacobian!(J, state_jacobian_type(TM), inplaceness(TM),
        model; kwargs...)
end

function get_state_jacobian!(J, model::TM, p; kwargs...) where TM
    return _get_state_jacobian!(J, state_jacobian_type(TM), inplaceness(TM),
        model, p; kwargs...)
end

function get_state_jacobian!(J, model::TM, dx, x, y, p, t; kwargs...) where TM
    return _get_state_jacobian!(J, state_jacobian_type(TM), inplaceness(TM),
        model, dx, x, y, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_state_jacobian!(J, ::Any, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_state_jacobian(model, args...; kwargs...)
end

# return an empty matrix
function _get_state_jacobian!(J, ::Empty, ::InPlace, models, args...; kwargs...)
    return J
end

# return a zero matrix
function _get_state_jacobian!(J, ::Zeros, ::InPlace, models, args...; kwargs...)
    return J .= 0
end

# return the identity matrix
function _get_state_jacobian!(J, ::Identity, ::InPlace, models, args...; kwargs...)
    J .= 0
    for i = 1:number_of_states(models)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_state_jacobian!` without arguments
function _get_state_jacobian!(J, ::Invariant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_state_jacobian!(J, models; kwargs...)
end

# dispatch to `get_state_jacobian!` without arguments
function _get_state_jacobian!(J, ::Invariant, ::InPlace, model, p; kwargs...)

    return get_state_jacobian!(J, model; kwargs...)
end

# dispatch to `get_state_jacobian!` with only parameters as arguments
function _get_state_jacobian!(J, ::Constant, ::InPlace, model, dx, x, y, p, t; kwargs...)

    return get_state_jacobian!(J, model, p; kwargs...)
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace, model, dx, x, y, p, t)

    f = x -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, x))
end

# calculate state jacobian for a combination of models
function _get_state_jacobian!(J, ::Invariant, ::InPlace, models::NTuple{N,AbstractModel};
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state and input indices
    ix = state_indices(models)
    iy = input_indices(models)

    # get input jacobian
    get_coupling_state_jacobian!(dydx_cache, models)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, models[i])
        for j = 1:length(models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, models[i])
                mul!(ri_xj, ri_yi, yi_xj, 1, 1)
            else
                mul!(ri_xj, ri_yi, yi_xj)
            end
        end
    end

    return J
end

# calculate state jacobian for a combination of models
function _get_state_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input state jacobian
    get_coupling_state_jacobian!(dydx_cache, models, p)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], pi)
        for j = 1:length(models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, models[i], pi)
                mul!(ri_xj, ri_yi, yi_xj, 1, 1)
            else
                mul!(ri_xj, ri_yi, yi_xj)
            end
        end
    end

    return J
end

# calculate state jacobian for a combination of models
function _get_state_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydx_cache = similar(J, number_of_inputs(models), number_of_states(models))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input state jacobian
    get_coupling_state_jacobian!(dydx_cache, models, dx, x, p, t)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], dxi, xi, yi, pi, t)
        for j = 1:length(models)
            ri_xj = view(J, ix[i], ix[j])
            yi_xj = view(dydx_cache, iy[i], ix[j])
            if i == j
                get_state_jacobian!(ri_xj, models[i], dxi, xi, yi, pi, t)
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
function get_input_jacobian(models::TM) where TM
    return _get_input_jacobian(input_jacobian_type(TM), inplaceness(TM), models)
end

function get_input_jacobian(models::TM, p) where TM
    return _get_input_jacobian(input_jacobian_type(TM), inplaceness(TM), models, p)
end

function get_input_jacobian(models::TM, dx, x, y, p, t) where TM
    return _get_input_jacobian(input_jacobian_type(TM), inplaceness(TM), models,
        dx, x, y, p, t)
end

# dispatch to an in-place function
function _get_input_jacobian(::Any, ::InPlace, models)
    Nx = number_of_states(models)
    Ny = number_of_inputs(models)
    J = Matrix{Float64}(undef, Nx, Ny)
    return get_input_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_input_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Nx = number_of_states(models)
    Ny = number_of_inputs(models)
    J = Matrix{TF}(undef, Nx, Ny)
    return get_input_jacobian!(J, models, p)
end

# dispatch to an in-place function
function _get_input_jacobian(::Any, ::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(models)
    Ny = number_of_inputs(models)
    J = Matrix{TF}(undef, Nx, Ny)
    return get_input_jacobian!(J, models, dx, x, y, p, t)
end

# return an empty matrix
function _get_input_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)
    return zero(SMatrix{Nx, Ny, Float64})
end

# return a zero matrix
function _get_input_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)
    return zero(SMatrix{Nx, Ny, Float64})
end

# return the identity matrix
function _get_input_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Ny = number_of_inputs(TM)
    return SMatrix{Nx, Ny, Float64}(I)
end

# dispatch to `get_input_jacobian` without arguments
function _get_input_jacobian(::Invariant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_input_jacobian(models)
end

# dispatch to the `get_input_jacobian` without arguments
function _get_input_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_input_jacobian(models)
end

# dispatch to `get_input_jacobian` with only parameters as arguments
function _get_input_jacobian(::Constant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_input_jacobian(models, p)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian(::Invariant, ::OutOfPlace, models)

    f = y -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, y)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian(::Constant, ::OutOfPlace, models, p)

    f = y -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, y)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace, models, dx, x, y, p, t)

    f = y -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, y)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    return get_invariant_input_partial(models)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, p) where N

    return get_constant_input_partial(models, p)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace, models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    return get_varying_input_partial(models, dx, x, y, p, t)
end

@generated function get_invariant_input_partial(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nx = number_of_states.(models.parameters)
    Ny = number_of_inputs.(models.parameters)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Ny[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_input_jacobian(models[1]), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_input_jacobian(models[$i]),
                $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_constant_input_partial(models::NTuple{N,AbstractModel},
    p) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Ny[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_input_jacobian(models[1], p[$(ip[1])]), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_input_jacobian(models[$i],
                p[$(ip[i])]), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

@generated function get_varying_input_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Ny[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_input_jacobian(models[1], dx[$(ix[1])], x[$(ix[1])],
            y[$(iy[1])], p[$(ip[1])], t), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_input_jacobian(models[$i],
                dx[$(ix[i])], x[$(ix[i])], y[$(iy[i])], p[$(ip[i])], t),
                $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

"""
    get_input_jacobian!(J, model)
    get_input_jacobian!(J, model, p)
    get_input_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_input_jacobian`](@ref)
"""
function get_input_jacobian!(J, model::TM; kwargs...) where TM
    return _get_input_jacobian!(J, input_jacobian_type(TM), inplaceness(TM),
        model; kwargs...)
end

function get_input_jacobian!(J, model::TM, p; kwargs...) where TM
    return _get_input_jacobian!(J, input_jacobian_type(TM), inplaceness(TM),
        model, p; kwargs...)
end

function get_input_jacobian!(J, model::TM, dx, x, y, p, t; kwargs...) where TM
    return _get_input_jacobian!(J, input_jacobian_type(TM), inplaceness(TM),
        model, dx, x, y, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_input_jacobian!(J, ::Any, ::OutOfPlace, model, args...; kwargs...)
    return J .= get_input_jacobian(model, args...; kwargs...)
end

# return an empty matrix
function _get_input_jacobian!(J, ::Empty, ::InPlace, model, args...; kwargs...)
    return J
end

# return a zero matrix
function _get_input_jacobian!(J, ::Zeros, ::InPlace, model, args...; kwargs...)
    return J .= 0
end

# return the identity matrix
function _get_input_jacobian!(J, ::Identity, ::InPlace, model, args...; kwargs...)
    J .= 0
    for i = 1:number_of_states(model)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_input_jacobian!` without arguments
function _get_input_jacobian!(J, ::Invariant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_input_jacobian!(J, models; kwargs...)
end

# dispatch to `get_input_jacobian!` without arguments
function _get_input_jacobian!(J, ::Invariant, ::InPlace, model, p;
    kwargs...)

    return get_input_jacobian!(J, models; kwargs...)
end

# dispatch to `get_input_jacobian!` with only parameters as arguments
function _get_input_jacobian!(J, ::Constant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_input_jacobian!(J, models, p; kwargs...)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian!(J, ::Invariant, ::InPlace, model)

    f = y -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, y))
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian!(J, ::Constant, ::InPlace, model, p)

    f = y -> get_residual(models, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, y))
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace, model, dx, x, y, p, t)

    f = y -> get_residual(model, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, y))
end

# calculate input jacobian for a combination of models
function _get_input_jacobian!(J, ::Invariant, ::InPlace, models::NTuple{N,AbstractModel}) where N

    # get state and input indices
    ix = state_indices(models)
    iy = input_indices(models)

    # calculate jacobian
    J .= 0
    for i = 1:length(models)
        Jii = view(J, ix[i], iy[j])
        get_state_jacobian!(Jii, models[i])
    end

    return J
end

# calculate input jacobian for a combination of models
function _get_input_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate jacobian
    J .= 0
    for i = 1:length(models)
        Jii = view(J, ix[i], iy[j])
        pi = view(p, ip[i])
        get_input_jacobian!(Jii, models[i], pi)
    end

    return J
end

# calculate input jacobian for a combination of models
function _get_input_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, y, p, t) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate jacobian
    J .= 0
    for i = 1:length(models)
        Jii = view(J, ix[i], iy[j])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        get_input_jacobian!(Jii, models[i], dxi, xi, yi, pi, t)
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
function get_parameter_jacobian(models::TM) where TM
    return _get_parameter_jacobian(parameter_jacobian_type(TM), inplaceness(TM), models)
end

function get_parameter_jacobian(models::TM, p) where TM
    return _get_parameter_jacobian(parameter_jacobian_type(TM), inplaceness(TM),
        models, p)
end

function get_parameter_jacobian(models::TM, dx, x, y, p, t) where TM
    return _get_parameter_jacobian(parameter_jacobian_type(TM), inplaceness(TM), models,
        dx, x, y, p, t)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::Any, ::InPlace, models)
    Nx = number_of_states(models)
    Np = number_of_parameters(models)
    J = Matrix{Float64}(undef, Nx, Np)
    return get_parameter_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Nx = number_of_states(models)
    Np = number_of_parameters(models)
    J = Matrix{TF}(undef, Nx, Np)
    return get_parameter_jacobian!(J, models, p)
end

# dispatch to an in-place function
function _get_parameter_jacobian(::Any, ::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(models)
    Np = number_of_parameters(models)
    J = Matrix{TF}(undef, Nx, Np)
    return get_parameter_jacobian!(J, models, dx, x, y, p, t)
end

# return an empty matrix
function _get_parameter_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Np = number_of_parameters(TM)
    return zero(SMatrix{Nx, Np, Float64})
end

# return a zero matrix
function _get_parameter_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Np = number_of_parameters(TM)
    return zero(SMatrix{Nx, Np, Float64})
end

# return the identity matrix
function _get_parameter_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    Np = number_of_parameters(TM)
    return SMatrix{Nx, Np, Float64}(I)
end

# dispatch to `get_parameter_jacobian` without arguments
function _get_parameter_jacobian(::Invariant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_parameter_jacobian(models)
end

# dispatch to `get_parameter_jacobian` without arguments
function _get_parameter_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_parameter_jacobian(models)
end

# dispatch to `get_parameter_jacobian` with only parameters as arguments
function _get_parameter_jacobian(::Constant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_parameter_jacobian(models, p)
end

# use automatic differentiation since a custom definition is absent
function _get_parameter_jacobian(::Invariant, ::OutOfPlace, models)

    f = p -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, p)
end

# use automatic differentiation since a custom definition is absent
function _get_parameter_jacobian(::Constant, ::OutOfPlace, models, p)

    f = p -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, p)
end

# use automatic differentiation since a custom definition is absent
function _get_parameter_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace, models, dx, x, y, p, t)

    f = p -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.jacobian(f, p)
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    r_p = get_parameter_partial(models)

    # calculate input jacobian
    r_y = get_input_jacobian(models)

    # calculate input jacobian
    y_p = get_coupling_parameter_jacobian(models)

    return r_p + r_y*y_p
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, p) where N

    # initialize parameter jacobian
    r_p = get_parameter_partial(models, p)

    # calculate input jacobian
    r_y = get_input_jacobian(models, p)

    # calculate coupling function parameter jacobian
    y_p = get_coupling_parameter_jacobian(models, p)

    return r_p + r_y*y_p
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, dx, x, y, p, t) where N

    # initialize parameter jacobian
    r_p = get_parameter_partial(models, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(models, dx, x, y, p, t)

    # calculate coupling function parameter jacobian
    y_p = get_coupling_parameter_jacobian(models, dx, x, p, t)

    return r_p + r_y*y_p
end

function get_parameter_partial(models::NTuple{N,AbstractModel}) where N
    get_invariant_parameter_partial(models)
end

function get_parameter_partial(models::NTuple{N,AbstractModel}, p) where N
    get_constant_parameter_partial(models, p)
end

function get_parameter_partial(models::NTuple{N,AbstractModel}, dx, x, y, p, t) where N
    get_varying_parameter_partial(models, dx, x, y, p, t)
end

@generated function get_invariant_parameter_partial(models::NTuple{N,AbstractModel}) where N

    # get number of state variables and parameters in each model
    Nx = number_of_states.(models.parameters)
    Np = number_of_parameters.(models.parameters)
    Npadd = number_of_additional_parameters(models.parameters...)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Np[j], Float64}) for i = 1:N, j = 1:N]

    # initialize columns for additional parameters
    Jadd = zeros(SMatrix{sum(Nx), Npadd})

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_parameter_jacobian(models[1]), $(Jij[1, 2:end]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[i, 1:i-1]...), get_parameter_jacobian(models[$i]), $(Jij[i, i+1:end]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    # add columns for additional parameters
    expr = quote
        $expr
        J = hcat(J, $Jadd)
    end

    return expr
end

@generated function get_constant_parameter_partial(models::NTuple{N,AbstractModel},
    p) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Np[j], Float64}) for i = 1:N, j = 1:N]

    # initialize additional columns for additional parameters
    Jadd = zeros(SMatrix{sum(Nx), Npadd})

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_parameter_jacobian(models[1], p[$(ip[1])]), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_parameter_jacobian(models[$i],
                p[$(ip[i])]), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    # add columns for additional parameters
    expr = quote
        $expr
        J = hcat(J, $Jadd)
    end

    return expr
end

@generated function get_varying_parameter_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nx[i], Np[j], Float64}) for i = 1:N, j = 1:N]

    # initialize additional columns for additional parameters
    Jadd = zeros(SMatrix{sum(Nx), Npadd})

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_parameter_jacobian(models[1], dx[$(ix[1])], x[$(ix[1])],
            y[$(iy[1])], p[$(ip[1])], t), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_parameter_jacobian(models[$i],
                dx[$(ix[i])], x[$(ix[i])], y[$(iy[i])], p[$(ip[i])], t), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    # add columns for additional parameters
    expr = quote
        $expr
        J = hcat(J, $Jadd)
    end

    return expr
end

"""
    get_parameter_jacobian!(J, model)
    get_parameter_jacobian!(J, model, p)
    get_parameter_jacobian!(J, model, dx, x, y, p, t)

In-place version of [`get_parameter_jacobian`](@ref)
"""
function get_parameter_jacobian!(J, model::TM; kwargs...) where TM
    return _get_parameter_jacobian!(J, parameter_jacobian_type(TM), inplaceness(TM),
        model; kwargs...)
end

function get_parameter_jacobian!(J, model::TM, p; kwargs...) where TM
    return _get_parameter_jacobian!(J, parameter_jacobian_type(TM), inplaceness(TM),
        model, p; kwargs...)
end

function get_parameter_jacobian!(J, model::TM, dx, x, y, p, t; kwargs...) where TM
    return _get_parameter_jacobian!(J, parameter_jacobian_type(TM), inplaceness(TM),
        model, dx, x, y, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_parameter_jacobian!(J, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return J .= get_parameter_jacobian(models, args...; kwargs...)
end

# return an empty matrix
function _get_parameter_jacobian!(J, ::Empty, ::InPlace, models, args...; kwargs...)
    return J
end

# return a zero matrix
function _get_parameter_jacobian!(J, ::Zeros, ::InPlace, models, args...; kwargs...)
    return J .= 0
end

# return the identity matrix
function _get_parameter_jacobian!(J, ::Identity, ::InPlace, models, args...; kwargs...)
    J .= 0
    for i = 1:number_of_states(models)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_parameter_jacobian!` without arguments
function _get_parameter_jacobian!(J, ::Invariant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_parameter_jacobian!(J, models; kwargs...)
end

# dispatch to `get_parameter_jacobian!` without arguments
function _get_parameter_jacobian!(J, ::Invariant, ::InPlace, model, p; kwargs...)

    return get_parameter_jacobian!(J, models; kwargs...)
end

# dispatch to `get_parameter_jacobian!` with only parameters as arguments
function _get_parameter_jacobian!(J, ::Constant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_parameter_jacobian!(J, models, p; kwargs...)
end

# use automatic differentiation since a custom definition is absent
function _get_parameter_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace, model, dx, x, y, p, t)

    f = p -> get_residual(model, dx, x, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, p))
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian!(J, ::Invariant, ::InPlace, models::NTuple{N,AbstractModel};
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydp_cache = similar(J, number_of_inputs(models), number_of_parameters(models))
    ) where N

    # get state and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input parameter jacobian
    get_coupling_parameter_jacobian!(dydp_cache, models)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, models[i])
        for j = 1:length(models)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, models[i])
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
    models::NTuple{N,AbstractModel}, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydp_cache = similar(J, number_of_inputs(models), number_of_parameters(models))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input jacobian
    get_coupling_parameter_jacobian!(dydp_cache, models, p)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], pi)
        for j = 1:length(models)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, models[i], pi)
                mul!(ri_pj, ri_yi, yi_pj, 1, 1)
            else
                mul!(ri_pj, ri_yi, yi_pj)
            end
        end
    end

    return J
end

# calculate parameter jacobian for a combination of models
function _get_parameter_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, y, p, t;
    drdy_cache = similar(J, number_of_states(models), number_of_inputs(models)),
    dydp_cache = similar(J, number_of_inputs(models), number_of_parameters(models))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input jacobian
    get_coupling_parameter_jacobian!(dydp_cache, models, dx, x, p, t)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], dxi, xi, yi, pi, t)
        for j = 1:length(models)
            ri_pj = view(J, ix[i], ip[j])
            yi_pj = view(dydp_cache, iy[i], ip[j])
            if i == j
                get_parameter_jacobian!(ri_pj, models[i], dxi, xi, yi, pi, t)
                mul!(ri_pj, ri_yi, yi_pj, 1, 1)
            else
                mul!(ri_pj, ri_yi, yi_pj)
            end
        end
    end

    return J
end

"""
    get_time_gradient(model, [dx, x, y, p, t])

Calculate the derivative of the residual function for `model` with respect to time
"""
function get_time_gradient(models::TM) where TM
    return _get_time_gradient(time_gradient_type(TM), inplaceness(TM), models)
end

function get_time_gradient(models::TM, p) where TM
    return _get_time_gradient(time_gradient_type(TM), inplaceness(TM), models, p)
end

function get_time_gradient(models::TM, dx, x, y, p, t) where TM
    return _get_time_gradient(time_gradient_type(TM), inplaceness(TM), models,
        dx, x, y, p, t)
end

# dispatch to an in-place function
function _get_time_gradient(::Any, ::InPlace, models)
    Nx = number_of_states(models)
    dT = Vector{Float64}(undef, Nx)
    return get_time_gradient!(dT, models)
end

# dispatch to an in-place function
function _get_time_gradient(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Nx = number_of_states(models)
    dT = Vector{TF}(undef, Nx)
    return get_time_gradient!(dT, models, p)
end

# dispatch to an in-place function
function _get_time_gradient(::Any, ::InPlace, models, dx, x, y, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(y), eltype(p), typeof(t))
    Nx = number_of_states(models)
    dT = Vector{TF}(undef, Nx)
    return get_time_gradient!(dT, models, dx, x, y, p, t)
end

# return an empty vector
function _get_time_gradient(::Empty, ::OutOfPlace, models, args...)
    return zeros(SVector{0, Float64})
end

# return a zero vector
function _get_time_gradient(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nx = number_of_states(TM)
    return zeros(SVector{Nx, Float64})
end

# dispatch to `get_time_gradient` without arguments
function _get_time_gradient(::Invariant, ::OutOfPlace, models, p)
    return get_time_gradient(models)
end

# dispatch to `get_time_gradient` without arguments
function _get_time_gradient(::Invariant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_time_gradient(models)
end

# dispatch to `get_time_gradient` with only parameters as arguments
function _get_time_gradient(::Constant, ::OutOfPlace, models, dx, x, y, p, t)
    return get_time_gradient(models, p)
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient(::Invariant, ::OutOfPlace, models)

    f = t -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.derivative(f, t)
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient(::Constant, ::OutOfPlace, models, p)

    f = t -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.derivative(f, t)
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient(::Union{Linear,Nonlinear}, ::OutOfPlace, models, dx, x, y, p, t)

    f = t -> get_residual(models, dx, x, y, p, t)

    return ForwardDiff.derivative(f, t)
end

# calculate time derivative for a combination of models
function _get_time_gradient(::Invariant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    r_t = get_time_partial(models)

    # calculate input jacobian
    r_y = get_input_jacobian(models)

    # calculate input jacobian
    y_t = get_coupling_time_gradient(models)

    return r_t + r_y*y_t
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient(::Constant, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, p) where N

    # initialize parameter jacobian
    r_t = get_time_partial(models, p)

    # calculate input jacobian
    r_y = get_input_jacobian(models, p)

    # calculate coupling function parameter jacobian
    y_t = get_coupling_time_gradient(models, p)

    return r_t + r_y*y_t
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, dx, x, y, p, t) where N

    # initialize parameter jacobian
    r_t = get_time_partial(models, dx, x, y, p, t)

    # calculate input jacobian
    r_y = get_input_jacobian(models, dx, x, y, p, t)

    # calculate coupling function parameter jacobian
    y_t = get_coupling_time_gradient(models, dx, x, p, t)

    return r_t + r_y*y_t
end

function get_time_partial(models::NTuple{N,AbstractModel}) where N
    return get_invariant_time_partial(models)
end

function get_time_partial(models::NTuple{N,AbstractModel}, p) where N
    return get_constant_time_partial(models, p)
end

function get_time_partial(models::NTuple{N,AbstractModel}, dx, x, y, p, t) where N
    return get_varying_time_partial(models, dx, x, y, p, t)
end

@generated function get_invariant_time_partial(models::NTuple{N,AbstractModel}) where N

    # initialize residual names
    dT = [Symbol("dT", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(dT[i]) = get_time_gradient(models[$i])
        end
    end

    expr = quote
        $expr
        dT = vcat($(dT...))
    end

    return expr
end

@generated function get_constant_time_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of parameters for each model
    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize residual names
    dT = [Symbol("dT", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(dT[i]) = get_time_gradient(models[$i], p[$(ip[i])])
        end
    end

    expr = quote
        $expr
        dT = vcat($(dT...))
    end

    return expr
end

@generated function get_varying_time_partial(models::NTuple{N,AbstractModel},
    dx, x, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nx = number_of_states.(models.parameters)
    ix2 = cumsum(Nx)
    ix1 = ix2 .- Nx .+ 1
    ix = ntuple(i->SVector{Nx[i]}(ix1[i]:ix2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize residual names
    dT = [Symbol("dT", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(dT[i]) = get_time_gradient(models[$i], dx[$(ix[i])], x[$(ix[i])],
                y[$(iy[i])], p[$(ip[i])], t)
        end
    end

    expr = quote
        $expr
        dT = vcat($(dT...))
    end

    return expr
end

"""
    get_time_gradient!(dT, model)
    get_time_gradient!(dT, model, p)
    get_time_gradient!(dT, model, dx, x, y, p, t)

In-place version of [`get_time_gradient`](@ref)
"""
function get_time_gradient!(dT, model::TM; kwargs...) where TM
    return _get_time_gradient!(dT, time_gradient_type(TM), inplaceness(TM),
        model; kwargs...)
end

function get_time_gradient!(dT, model::TM, p; kwargs...) where TM
    return _get_time_gradient!(dT, time_gradient_type(TM), inplaceness(TM),
        model, p; kwargs...)
end

function get_time_gradient!(dT, model::TM, dx, x, y, p, t; kwargs...) where TM
    return _get_time_gradient!(dT, time_gradient_type(TM), inplaceness(TM),
        model, dx, x, y, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_time_gradient!(dT, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return dT .= get_time_gradient(models, args...; kwargs...)
end

# return an empty vector
function _get_time_gradient!(dT, ::Empty, ::InPlace, models, args...; kwargs...)
    return dT
end

# return a zero vector
function _get_time_gradient!(dT, ::Zeros, ::InPlace, models, args...; kwargs...)
    return dT .= 0
end

# dispatch to `get_time_gradient!` without arguments
function _get_time_gradient!(dT, ::Invariant, ::InPlace, model, p;
    kwargs...)

    return get_time_gradient!(dT, models; kwargs...)
end

# dispatch to `get_time_gradient!` without arguments
function _get_time_gradient!(dT, ::Invariant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_time_gradient!(dT, models; kwargs...)
end

# dispatch to `get_time_gradient!` with only parameters as arguments
function _get_time_gradient!(dT, ::Constant, ::InPlace, model, dx, x, y, p, t;
    kwargs...)

    return get_time_gradient!(dT, models, p; kwargs...)
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient!(dT, ::Invariant, ::InPlace, model)

    f = t -> get_residual(models, dx, x, y, p, t)

    return copyto!(dT, ForwardDiff.derivative(f, t))
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient!(dT, ::Constant, ::InPlace, model, p)

    f = t -> get_residual(models, dx, x, y, p, t)

    return copyto!(dT, ForwardDiff.derivative(f, t))
end

# use automatic differentiation since a custom definition is absent
function _get_time_gradient!(dT, ::Union{Linear,Nonlinear}, ::InPlace, model, dx, x, y, p, t)

    f = t -> get_residual(models, dx, x, y, p, t)

    return copyto!(dT, ForwardDiff.derivative(f, t))
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient!(dT, ::Invariant, ::InPlace, models::NTuple{N,AbstractModel};
    drdy_cache = similar(dT, number_of_states(models), number_of_inputs(models)),
    dydt_cache = similar(dT, number_of_inputs(models))
    ) where N

    # get state and input indices
    ix = state_indices(models)
    iy = input_indices(models)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, models)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        get_input_jacobian!(ri_yi, models[i])
        for j = 1:length(models)
            ri_t = view(dT, ix[i])
            yi_t = view(dydt_cache, iy[i])
            if i == j
                get_parameter_jacobian!(Jij, models[i])
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient!(dT, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p;
    drdy_cache = similar(dT, number_of_states(models), number_of_inputs(models)),
    dydt_cache = similar(dT, number_of_inputs(models))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, models, p)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], pi)
        for j = 1:length(models)
            ri_t = view(dT, ix[i], ix[j])
            yi_t = view(dydt_cache, iy[i], ix[j])
            if i == j
                get_parameter_jacobian!(ri_t, models[i], pi)
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

# calculate parameter jacobian for a combination of models
function _get_time_gradient!(dT, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, y, p, t;
    drdy_cache = similar(dT, number_of_states(models), number_of_inputs(models)),
    dydt_cache = similar(dT, number_of_inputs(models))
    ) where N

    # get state, input, and parameter indices
    ix = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # get input jacobian
    get_coupling_time_gradient!(dydt_cache, models, dx, x, p, t)

    # calculate jacobian
    for i = 1:length(models)
        ri_yi = view(drdy_cache, ix[i], iy[i])
        dxi = view(dx, ix[i])
        xi = view(x, ix[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        ri_yi = get_input_jacobian!(ri_yi, models[i], dxi, xi, yi, pi, t)
        for j = 1:length(models)
            ri_t = view(dT, ix[i], ix[j])
            yi_t = view(dydt_cache, iy[i], ix[j])
            if i == j
                get_parameter_jacobian!(ri_t, models[i], dxi, xi, yi, pi, t)
                mul!(ri_t, ri_yi, yi_t, 1, 1)
            else
                mul!(ri_t, ri_yi, yi_t)
            end
        end
    end

    return dT
end

"""
    get_coupling_inputs(models, dx, x, p, t)

Calculate the inputs to the specified combination of models.
"""
function get_coupling_inputs(models::TM, dx, x, p, t) where TM
    return _get_coupling_inputs(coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# dispatch to an in-place function
function _get_coupling_inputs(::InPlace, models, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(TF, Ny)
    get_coupling_inputs!(y, models, dx, x, p, t)
    return y
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_inputs(::OutOfPlace, models::NTuple{N,AbstractModel}, dx, x, p, t) where N
    return get_coupling_inputs(models..., dx, x, p, t)
end

"""
    get_coupling_inputs!(y, models::NTuple{N,AbstractModel}, dx, x, p, t) where N

In-place version of [`get_coupling_inputs`](@ref)
"""
function get_coupling_inputs!(y, models::T, dx, x, p, t) where T
    return _get_coupling_inputs!(y, coupling_inplaceness(T.parameters...), models,
        dx, x, p, t)
end

# dispatch to an out-of-place function
function _get_coupling_inputs!(y, ::OutOfPlace, models, dx, x, p, t)
    return y .= get_coupling_inputs(models, dx, x, p, t)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_inputs!(y, ::InPlace, models::NTuple{N,AbstractModel}, dx, x, p, t) where N
    return get_coupling_inputs!(y, models..., dx, x, p, t)
end

"""
    get_coupling_rate_jacobian(model)
    get_coupling_rate_jacobian(model, p)
    get_coupling_rate_jacobian(model, dx, x, p, t)

Calculate the jacobian of the coupling function for `model` with respect to the
state rates
"""
get_coupling_rate_jacobian

function get_coupling_rate_jacobian(models::TM) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian(coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_rate_jacobian(models::TM, p) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian(coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_rate_jacobian(models::TM, dx, x, p, t) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian(coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if the jacobian is not defined
function get_coupling_rate_jacobian(args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = dx -> get_coupling_inputs(models, dx, x, p, t)

    return ForwardDiff.jacobian(f, dx)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::Any, ::InPlace, models)
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{Float64}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, models, p)
end

# dispatch to an in-place function
function _get_coupling_rate_jacobian(::Any, ::InPlace, models, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_rate_jacobian!(J, models, dx, x, p, t)
end

# return an empty matrix
function _get_coupling_rate_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return zeros(SMatrix{Ny, Nx, Float64})
end

# return a zero matrix
function _get_coupling_rate_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return zeros(SMatrix{Ny, Nx, Float64})
end

# return the identity matrix
function _get_coupling_rate_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return SMatrix{Ny, Nx, Float64}(I)
end

# dispatch to the `get_coupling_rate_jacobian` function without arguments
function _get_coupling_rate_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_coupling_rate_jacobian(models)
end

# dispatch to the `get_coupling_rate_jacobian` function without arguments
function _get_coupling_rate_jacobian(::Invariant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_rate_jacobian(models)
end

# dispatch to `get_coupling_rate_jacobian` with only parameters as arguments
function _get_coupling_rate_jacobian(::Constant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_rate_jacobian(models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}) where N

    return get_coupling_rate_jacobian(models...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, p) where N

    return get_coupling_rate_jacobian(models..., p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    return get_coupling_rate_jacobian(models..., dx, x, p, t)
end

"""
    get_coupling_rate_jacobian!(J, model)
    get_coupling_rate_jacobian!(J, model, p)
    get_coupling_rate_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_rate_jacobian`](@ref).
"""
get_coupling_rate_jacobian!

function get_coupling_rate_jacobian!(J, models::TM) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian!(J, coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_rate_jacobian!(J, models::TM, p) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian!(J, coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_rate_jacobian!(J, models::TM, dx, x, p, t) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_rate_jacobian!(J, coupling_rate_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if the jacobian is not defined
function get_coupling_rate_jacobian!(J, args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = dx -> get_coupling_inputs(models, dx, x, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, dx))
end

# dispatch to an out-of-place function
function _get_coupling_rate_jacobian!(J, ::Any, ::OutOfPlace, models, args...)
    return J .= get_coupling_rate_jacobian(models, args...)
end

# return an empty matrix
function _get_coupling_rate_jacobian!(J, ::Empty, ::InPlace, models, args...)
    return J
end

# return a zero matrix
function _get_coupling_rate_jacobian!(J, ::Zeros, ::InPlace, models, args...)
    return J .= 0
end

# return the identity matrix
function _get_coupling_rate_jacobian!(J, ::Identity, ::InPlace, models, args...)
    J .= 0
    for i = 1:number_of_inputs(models)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_coupling_rate_jacobian!` without arguments
function _get_coupling_rate_jacobian!(J, ::Invariant, ::InPlace, models, p)

    return get_coupling_rate_jacobian!(J, models)
end

# dispatch to `get_coupling_rate_jacobian!` without arguments
function _get_coupling_rate_jacobian!(J, ::Invariant, ::InPlace, models, dx, x, p, t)

    return get_coupling_rate_jacobian!(J, models)
end

# dispatch to `get_coupling_rate_jacobian!` with only parameters as arguments
function _get_coupling_rate_jacobian!(J, ::Constant, ::InPlace, models, dx, x, p, t)

    return get_coupling_rate_jacobian!(J, models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian!(J, ::Invariant, ::InPlace,
    models::NTuple{N,AbstractModel}) where N

    get_coupling_rate_jacobian!(J, models...)

    return J
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p) where N

    get_coupling_rate_jacobian!(J, models..., p)

    return J
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_rate_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    get_coupling_rate_jacobian!(J, models..., dx, x, p, t)

    return J
end

"""
    get_coupling_state_jacobian(model)
    get_coupling_state_jacobian(model, p)
    get_coupling_state_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the coupling function for `model` with respect to the
state variables
"""
get_coupling_state_jacobian

function get_coupling_state_jacobian(models::TM) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_state_jacobian(coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_state_jacobian(models::TM, p) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_state_jacobian(coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_state_jacobian(models::TM, dx, x, p, t) where TM <: NTuple{N, AbstractModel} where N
    return _get_coupling_state_jacobian(coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_state_jacobian(args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = x -> get_coupling_inputs(models, dx, x, p, t)

    return ForwardDiff.jacobian(f, x)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::Any, ::InPlace, models)
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{Float64}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, models, dx, x, p, t)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::Any, ::InPlace, models, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    Nx = number_of_states(models)
    J = Matrix{TF}(undef, Ny, Nx)
    return get_coupling_state_jacobian!(J, models, dx, x, p, t)
end

# return an empty matrix
function _get_coupling_state_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return zeros(SMatrix{Ny, Nx, Float64})
end

# return a matrix of zeros
function _get_coupling_state_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return zeros(SMatrix{Ny, Nx, Float64})
end

# return the identity matrix
function _get_coupling_state_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Nx = number_of_states(TM)
    return SMatrix{Ny, Nx, Float64}(I)
end

# dispatch to `get_coupling_state_jacobian` function without arguments
function _get_coupling_state_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_coupling_state_jacobian(models)
end

# dispatch to `get_coupling_state_jacobian` function without arguments
function _get_coupling_state_jacobian(::Invariant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_state_jacobian(models)
end

# dispatch to `get_coupling_state_jacobian` with only parameters as arguments
function _get_coupling_state_jacobian(::Constant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_state_jacobian(models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}) where N

    return get_coupling_state_jacobian(models...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, p) where N

    return get_coupling_state_jacobian(models..., p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    return get_coupling_state_jacobian(models..., dx, x, p, t)
end

"""
    get_coupling_state_jacobian!(J, model)
    get_coupling_state_jacobian!(J, model, p)
    get_coupling_state_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_state_jacobian`](@ref)
"""
get_coupling_state_jacobian!

function get_coupling_state_jacobian!(J, models::TM) where TM
    return _get_coupling_state_jacobian!(J, coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_state_jacobian!(J, models::TM, p) where TM
    return _get_coupling_state_jacobian!(J, coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_state_jacobian!(J, models::TM, dx, x, p, t) where TM
    return _get_coupling_state_jacobian!(J, coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_state_jacobian!(J, args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = x -> get_coupling_inputs(models, dx, x, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, x))
end

# dispatch to an out-of-place function
function _get_coupling_state_jacobian!(J, ::Any, ::OutOfPlace, models, args...)
    return J .= get_coupling_state_jacobian(models, args...)
end

# return an empty matrix
function _get_coupling_state_jacobian!(J, ::Empty, ::InPlace, models, args...)
    return J
end

# return a zero matrix
function _get_coupling_state_jacobian!(J, ::Zeros, ::InPlace, models, args...)
    return J .= 0
end

# return the identity matrix
function _get_coupling_state_jacobian!(J, ::Identity, ::InPlace, models, args...)
    J .= 0
    for i = 1:number_of_inputs(models)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_coupling_state_jacobian!` without arguments
function _get_coupling_state_jacobian!(J, ::Invariant, ::InPlace, models, p)

    return get_coupling_state_jacobian!(J, models)
end

# dispatch to `get_coupling_state_jacobian!` without arguments
function _get_coupling_state_jacobian!(J, ::Invariant, ::InPlace, models, dx, x, p, t)

    return get_coupling_state_jacobian!(J, models)
end

# dispatch to `get_coupling_rate_jacobian!` with only parameters as arguments
function _get_coupling_state_jacobian!(J, ::Constant, ::InPlace, models, dx, x, p, t)

    return get_coupling_state_jacobian!(J, models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian!(J, ::Invariant, ::InPlace,
    models::NTuple{N,AbstractModel}) where N

    get_coupling_state_jacobian!(J, models...)

    return J
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p) where N

    get_coupling_state_jacobian!(J, models..., p)

    return J
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    get_coupling_state_jacobian!(J, models..., dx, x, p, t)

    return J
end

"""
    get_coupling_parameter_jacobian(model)
    get_coupling_parameter_jacobian(model, p)
    get_coupling_parameter_jacobian(model, dx, x, y, p, t)

Calculate the jacobian of the coupling function for `model` with respect to the
parameters
"""
get_coupling_parameter_jacobian

function get_coupling_parameter_jacobian(models::TM) where TM
    return _get_coupling_parameter_jacobian(coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_parameter_jacobian(models::TM, p) where TM
    return _get_coupling_parameter_jacobian(coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_parameter_jacobian(models::TM, dx, x, p, t) where TM
    return _get_coupling_parameter_jacobian(coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_parameter_jacobian(args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = p -> get_coupling_inputs(models, dx, x, p, t)

    return ForwardDiff.jacobian(f, p)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::Any, ::InPlace, models)
    Ny = number_of_inputs(models)
    Np = number_of_parameters(models)
    J = Matrix{Float64}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Ny = number_of_inputs(models)
    Np = number_of_parameters(models)
    J = Matrix{TF}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, models)
end

# dispatch to an in-place function
function _get_coupling_parameter_jacobian(::Any, ::InPlace, models, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    Np = number_of_parameters(models)
    J = Matrix{TF}(undef, Ny, Np)
    return get_coupling_parameter_jacobian!(J, models, dx, x, p, t)
end

# return an empty matrix
function _get_coupling_parameter_jacobian(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Np = number_of_parameters(TM)
    return zeros(SMatrix{Ny, Np, Float64})
end

# return a zero matrix
function _get_coupling_parameter_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Np = number_of_parameters(TM)
    return zeros(SMatrix{Ny, Np, Float64})
end

# return the identity matrix
function _get_coupling_parameter_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    Np = number_of_parameters(TM)
    return SMatrix{Ny, Np, Float64}(I)
end

# dispatch to the `get_coupling_parameter_jacobian` function without arguments
function _get_coupling_parameter_jacobian(::Invariant, ::OutOfPlace, models, p)
    return get_coupling_parameter_jacobian(models)
end

# dispatch to the `get_coupling_parameter_jacobian` function without arguments
function _get_coupling_parameter_jacobian(::Invariant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_parameter_jacobian(models)
end

# dispatch to `get_coupling_parameter_jacobian` with only parameters as arguments
function _get_coupling_parameter_jacobian(::Constant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_parameter_jacobian(models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian(::Invariant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}) where N

    return get_coupling_parameter_jacobian(models...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, p) where N

    return get_coupling_parameter_jacobian(models..., p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    return get_coupling_parameter_jacobian(models..., dx, x, p, t)
end

"""
    get_coupling_parameter_jacobian!(J, model)
    get_coupling_parameter_jacobian!(J, model, p)
    get_coupling_parameter_jacobian!(J, model, dx, x, p, t)

In-place version of [`get_coupling_parameter_jacobian`](@ref)
"""
get_coupling_parameter_jacobian!

function get_coupling_parameter_jacobian!(J, models::TM) where TM
    return _get_coupling_parameter_jacobian!(J, coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_parameter_jacobian!(J, models::TM, p) where TM
    return _get_coupling_parameter_jacobian!(J, coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_parameter_jacobian!(J, models::TM, dx, x, p, t) where TM
    return _get_coupling_parameter_jacobian!(J, coupling_parameter_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_parameter_jacobian!(J, args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = p -> get_coupling_inputs(models, dx, x, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, p))
end

# dispatch to an out-of-place function
function _get_coupling_parameter_jacobian!(J, ::Any, ::OutOfPlace, models, args...)
    return J .= get_coupling_parameter_jacobian(models, args...)
end

# return an empty matrix
function _get_coupling_parameter_jacobian!(J, ::Empty, ::InPlace, models, args...)
    return J
end

# return a zero matrix
function _get_coupling_parameter_jacobian!(J, ::Zeros, ::InPlace, models, args...)
    return J .= 0
end

# return the identity matrix
function _get_coupling_parameter_jacobian!(J, ::Identity, ::InPlace, models, args...)
    J .= 0
    for i = 1:number_of_inputs(models)
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_coupling_parameter_jacobian!` without arguments
function _get_coupling_parameter_jacobian!(J, ::Invariant, ::InPlace, models, p)
    return get_coupling_parameter_jacobian!(J, models)
end

# dispatch to `get_coupling_parameter_jacobian!` without arguments
function _get_coupling_parameter_jacobian!(J, ::Invariant, ::InPlace, models, dx, x, p, t)
    return get_coupling_parameter_jacobian!(J, models)
end

# dispatch to `get_coupling_parameter_jacobian!` with only parameters as arguments
function _get_coupling_parameter_jacobian!(J, ::Constant, ::InPlace, models, dx, x, p, t)
    return get_coupling_parameter_jacobian!(J, models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian!(J, ::Invariant, ::InPlace,
    models::NTuple{N,AbstractModel}) where N

    get_coupling_parameter_jacobian!(J, models...)

    return J
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian!(J, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p) where N

    get_coupling_parameter_jacobian!(J, models..., p)

    return J
end


# dispatch to the user-provided function for the specific combination of models
function _get_coupling_parameter_jacobian!(J, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    get_coupling_parameter_jacobian!(J, models..., dx, x, p, t)

    return J
end

"""
    get_coupling_time_gradient(model)
    get_coupling_time_gradient(model, p)
    get_coupling_time_gradient(model, dx, x, y, p, t)

Calculate the derivative of the coupling function for `model` with respect to time
"""
get_coupling_time_gradient

function get_coupling_time_gradient(models::TM) where TM
    return _get_coupling_time_gradient(coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_time_gradient(models::TM, p) where TM
    return _get_coupling_time_gradient(coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_time_gradient(models::TM, dx, x, p, t) where TM
    return _get_coupling_time_gradient(coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_time_gradient(args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = t -> get_coupling_inputs(models, dx, x, p, t)

    return ForwardDiff.derivative(f, t)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::Any, ::InPlace, models)
    Ny = number_of_inputs(models)
    dT = Vector{Float64}(undef, Ny)
    return get_coupling_time_gradient!(dT, models)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::Any, ::InPlace, models, p)
    TF = eltype(p)
    Ny = number_of_inputs(models)
    dT = Vector{TF}(undef, Ny)
    return get_coupling_time_gradient!(dT, models, p)
end

# dispatch to an in-place function
function _get_coupling_time_gradient(::Any, ::InPlace, models, dx, x, p, t)
    TF = promote_type(eltype(dx), eltype(x), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    dT = Vector{TF}(undef, Ny)
    return get_coupling_time_gradient!(dT, models, dx, x, p, t)
end

# return an empty vector
function _get_coupling_time_gradient(::Empty, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    return zeros(SVector{Ny, Float64})
end

# return a zero vector
function _get_coupling_time_gradient(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    return zeros(SVector{Ny, Float64})
end

# dispatch to the `get_coupling_time_gradient` function without arguments
function _get_coupling_time_gradient(::Invariant, ::OutOfPlace, models, p)
    return get_coupling_time_gradient(models)
end

# dispatch to the `get_coupling_time_gradient` function without arguments
function _get_coupling_time_gradient(::Invariant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_time_gradient(models)
end

# dispatch to `get_coupling_time_gradient` with only parameters as arguments
function _get_coupling_time_gradient(::Constant, ::OutOfPlace, models, dx, x, p, t)
    return get_coupling_time_gradient(models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient(::Invariant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}) where N

    return get_coupling_time_gradient(models...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, p) where N

    return get_coupling_time_gradient(models..., p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient(::Union{Linear,Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    return get_coupling_time_gradient(models..., dx, x, p, t)
end

"""
    get_coupling_time_gradient!(dT, model)
    get_coupling_time_gradient!(dT, model, p)
    get_coupling_time_gradient!(dT, model, dx, x, p, t)

In-place version of [`get_coupling_time_gradient`](@ref)
"""
get_coupling_time_gradient!

function get_coupling_time_gradient!(dT, models::TM) where TM
    return _get_coupling_time_gradient!(dT, coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models)
end

function get_coupling_time_gradient!(dT, models::TM, p) where TM
    return _get_coupling_time_gradient!(dT, coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, p)
end

function get_coupling_time_gradient!(dT, models::TM, dx, x, p, t) where TM
    return _get_coupling_time_gradient!(dT, coupling_time_gradient_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, dx, x, p, t)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_time_gradient!(dT, args...)

    models = args[1:end-4]
    dx = args[end-3]
    x = args[end-2]
    p = args[end-1]
    t = args[end]

    f = t -> get_coupling_inputs(models, dx, x, p, t)

    return copyto!(dT, ForwardDiff.derivative(f, t))
end

# dispatch to an out-of-place function
function _get_coupling_time_gradient!(dT, ::Any, ::OutOfPlace, models, args...)
    return dT .= get_coupling_time_gradient(models, args...)
end

# return an empty matrix
function _get_coupling_time_gradient!(dT, ::Empty, ::InPlace, models, args...)
    return dT
end

# return a matrix of zeros
function _get_coupling_time_gradient!(dT, ::Zeros, ::InPlace, models, args...)
    return dT .= 0
end

# dispatch to `get_coupling_time_gradient!` without arguments
function _get_coupling_time_gradient!(dT, ::Invariant, ::InPlace, models, p)
    return get_coupling_time_gradient!(dT, models)
end

# dispatch to `get_coupling_time_gradient!` without arguments
function _get_coupling_time_gradient!(dT, ::Invariant, ::InPlace, models, dx, x, p, t)
    return get_coupling_time_gradient!(dT, models)
end

# dispatch to `get_coupling_time_gradient!` with only parameters as arguments
function _get_coupling_time_gradient!(dT, ::Constant, ::InPlace, models, dx, x, p, t)
    return get_coupling_time_gradient!(dT, models, p)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient!(dT, ::Invariant, ::InPlace,
    models::NTuple{N,AbstractModel}) where N

    get_coupling_time_gradient!(dT, models...)

    return dT
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient!(dT, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}, p) where N

    get_coupling_time_gradient!(dT, models..., p)

    return dT
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_time_gradient!(dT, ::Union{Linear,Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    get_coupling_time_gradient!(dT, models..., dx, x, p, t)

    return dT
end

"""
    get_eigen(model::TM, x, p; kwargs...)

Return the eigenvalues, left eigenvector matrix, and right eigenvector matrix
corresponding to the model for state variables `x` and parameters `p`.

For in-place models, the number of eigenvalues to compute may be specified using
the `nev` keyword argument.
"""
get_eigen

function get_eigen(model::TM, x, p; dx = FillArrays.Zeros(x), t = 0,
    kwargs...) where TM <: AbstractModel

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = 1 : Np
    iy = Np + 1 : Np + Ny

    return _get_eigen(inplaceness(TM), model, dx, x, view(p,ip), view(p,iy), t; kwargs...)
end

function get_eigen(model::TM, x, p; dx = FillArrays.Zeros(x), t = 0,
    y = get_coupling_inputs(model, dx, x, p, t), kwargs...) where TM <: NTuple{N,AbstractModel} where N

    return _get_eigen(inplaceness(TM), model, dx, x, y, p, t; kwargs...)
end

function _get_eigen(::OutOfPlace, model, args...)
    A = -Array(get_state_jacobian(model, args...)) # jacobian
    B = Array(get_rate_jacobian(model, args...)) # mass matrix
    E = eigen(A, B) # eigenvalue decomposition
    λ = eigvals(E) # eigenvalues
    V = eigvecs(E) # right eigenvector matrix
    U = I/(A*V) # left eigenvector matrix
    return λ, U, V
end

function _get_eigen(::InPlace, model, args...; nev=min(20, number_of_states(model)))

    # calculate state and rate jacobians
    A = -get_state_jacobian(model, args...)
    B = get_rate_jacobian(model, args...)

    # construct linear map
    nx = size(A, 1)
    TF = promote_type(eltype(A), eltype(B))
    Afact = lu(A)
    f! = (b, x) -> ldiv!(b, Afact, B*x)
    fc! = (b, x) -> mul!(b, B', Afact'\x)
    Abar = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

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

"""
    get_ode(model)
    get_ode(model, p)

Construct an ODEFunction corresponding to the specified model or models which
may be solved using DifferentialEquations.
"""
function get_ode(model::TM, args...) where TM
    return _get_ode(rate_jacobian_type(TM), inplaceness(TM), model, args...)
end

function _get_ode(::Identity, ::OutOfPlace, model::AbstractModel, p=nothing)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dx .*= -1
    end
    mass_matrix = I
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Identity, ::OutOfPlace, model::NTuple{N,AbstractModel}, p=nothing) where N

    fy = (x, p, t) -> get_coupling_inputs(models, FillArrays.Zeros(x), x, p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = I
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Invariant, ::OutOfPlace, model::AbstractModel, p=nothing)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(models)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, models, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Invariant, ::OutOfPlace, model::NTuple{N,AbstractModel}, p=nothing) where N

    fy = (x, p, t) -> get_coupling_inputs(models, FillArrays.Zeros(x), x, p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, models, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

_get_ode(::Constant, ::OutOfPlace, model) = _get_ode(Linear(), OutOfPlace(), model)

function _get_ode(::Constant, ::OutOfPlace, model::AbstractModel, p)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, x, p[iy], p[ip], t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model, p)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Constant, ::OutOfPlace, model::NTuple{N,AbstractModel}, p) where N

    fy = (x, p, t) -> get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model, p)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Linear, ::OutOfPlace, model::AbstractModel, p=nothing)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_rate_jacobian!(M, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Linear, ::OutOfPlace, model::NTuple{N,AbstractModel}, p=nothing) where N

    Nx = number_of_states(model)

    fy = (x, p, t) -> get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_rate_jacobian!(M, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Identity, ::InPlace, model::AbstractModel, p=nothing)

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dx .*= -1
    end
    mass_matrix = I
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Identity, ::InPlace, model::NTuple{N,AbstractModel}, p=nothing) where N

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # cached variables
    ucache = fill(NaN, Nx)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    fy = (x, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
            # update the cache variables
            if (x != ucache) && (p != pcache) && (t != tcache[])
                ucache .= x
                ycache .= get_coupling_inputs!(ycache, model, FillArrays.Zeros(x), x, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
        end
        return y
    end

    # TODO: use cached variables for input jacobian and coupling function jacobians

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = I
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Invariant, ::InPlace, model::AbstractModel, p=nothing)

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Invariant, ::InPlace, model::NTuple{N,AbstractModel}, p=nothing) where N

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # cached variables
    ucache = fill(NaN, Nx)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    fy = (x, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
            # update the cache variables
            if (x != ucache) && (p != pcache) && (t != tcache[])
                ucache .= x
                ycache .= get_coupling_inputs!(ycache, model, FillArrays.Zeros(x), x, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
        end
        return y
    end

    # TODO: use cached variables for input jacobian and coupling function jacobians

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Constant, ::InPlace, model::AbstractModel)

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_rate_jacobian!(M, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Constant, ::InPlace, model::NTuple{N,AbstractModel}) where N

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # cached variables
    ucache = fill(NaN, Nx)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    fy = (x, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
            # update the cache variables
            if (x != ucache) && (p != pcache) && (t != tcache[])
                ucache .= x
                ycache .= get_coupling_inputs!(ycache, model, FillArrays.Zeros(x), x, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
        end
        return y
    end

    # TODO: use cached variables for input jacobian and coupling function jacobians

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_rate_jacobian!(M, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end


function _get_ode(::Constant, ::InPlace, model::AbstractModel, p)

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model, p)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Constant, ::InPlace, model::NTuple{N,AbstractModel}, p) where N

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # cached variables
    ucache = fill(NaN, Nx)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    fy = (x, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
            # update the cache variables
            if (x != ucache) && (p != pcache) && (t != tcache[])
                ucache .= x
                ycache .= get_coupling_inputs!(ycache, model, FillArrays.Zeros(x), x, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
        end
        return y
    end

    # TODO: use cached variables for input jacobian and coupling function jacobians

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = get_rate_jacobian(model, p)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Linear, ::InPlace, model::AbstractModel, p=nothing)

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_rate_jacobian!(M, model, FillArrays.Zeros(x), x, p[iy], p[ip], t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, view(p, iy), view(p, ip), t)
        pJ .*= -1
    end
    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

function _get_ode(::Linear, ::InPlace, model::NTuple{N,AbstractModel}, p) where N

    # problem dimensions
    Nx = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # cached variables
    ucache = fill(NaN, Nx)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    fy = (x, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if promote_type(eltype(x), eltype(p), typeof(t)) <: Float64
            # update the cache variables
            if (x != ucache) && (p != pcache) && (t != tcache[])
                ucache .= x
                ycache .= get_coupling_inputs!(ycache, model, FillArrays.Zeros(x), x, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_coupling_inputs(model, FillArrays.Zeros(x), x, p, t)
        end
        return y
    end

    # TODO: use cached variables for input jacobian and coupling function jacobians

    M = zeros(Nx, Nx)
    update_func = (M, x, p, t) -> get_state_jacobian!(M, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dx .*= -1
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)
    tgrad = (dT, x, p, t) -> begin
        get_time_gradient!(dT, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        dT .*= -1
    end
    jac = (J, x, p, t) -> begin
        get_state_jacobian!(J, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        J .*= -1
    end
    paramjac = (pJ, x, p, t) -> begin
        get_parameter_jacobian!(pJ, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        pJ .*= -1
    end

    # construct and return an ODEFunction
    return ODEFunction{true,true}(f; mass_matrix, tgrad, jac, paramjac)
end

# plotting recipes
@recipe function f(models::AbstractModel, dx, x, p, t) where N

    Np = number_of_parameters(models)
    Ny = number_of_inputs(models)

    ip = 1:Np
    iy = Np+1:Np+Ny

    return models..., dx, x, p[iy], p[ip], t
end

@recipe function f(models::AbstractModel, x, p, t) where N

    return models, FillArrays.Zeros(x), x, p, t
end

@recipe function f(models::AbstractModel, sol, t)

    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1 : Np
    iy = Np + 1 : Np + Ny

    x = sol(t)
    y = sol.prob.p[iy]
    p = sol.prob.p[ip]

    dx = sol(t, Val{1})

    return models..., dx, x, y, p, t
end

@recipe function f(model::AbstractModel, sol)

    it = sol.tslocation

    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1 : Np
    iy = Np + 1 : Np + Ny

    x = sol.x[it]
    y = sol.prob.p[iy]
    p = sol.prob.p[ip]
    t = sol.t[it]

    dx = sol(t, Val{1})

    return model, dx, x, y, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, dx, x, p, t) where N

    y = get_coupling_inputs(models, dx, x, p, t)

    return models..., dx, x, y, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, x, p, t) where N

    return models, FillArrays.Zeros(x), x, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, sol, t) where N

    x = sol(t)
    p = sol.prob.p

    dx = sol(t, Val{1})

    y = get_coupling_inputs(models, dx, x, p, t)

    return models..., dx, x, y, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, sol) where N

    it = sol.tslocation

    x = sol.x[it]
    p = sol.prob.p
    t = sol.t[it]

    dx = sol(t, Val{1})

    y = get_coupling_inputs(models, dx, x, p, t)

    return models..., dx, x, y, p, t
end
