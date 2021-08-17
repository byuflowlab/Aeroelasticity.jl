"""
    couple_models(models...)

Function which couples multiple models together to form a coupled model.
"""
couple_models(models...) = models

"""
    number_of_states(model)

Return the total number of states corresponding to the model or models.
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
function number_of_states(::Type{TM}) where TM<:NTuple{N,AbstractModel} where N
    return sum(number_of_states.(TM.parameters))
end

"""
    number_of_inputs(model)

Return the total number of inputs corresponding to the model or models.
"""
number_of_inputs

# by default, dispatch on type
number_of_inputs(model::TM) where TM <: AbstractModel = number_of_inputs(TM)

# models with no states have no inputs
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

# coupled models have concatenated parameters
function number_of_parameters(models::NTuple{N,AbstractModel}) where N
    sum(number_of_parameters.(models)) + number_of_additional_parameters(models...)
end

# coupled models have concatenated parameters
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
    Nu = number_of_states.(models)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    return UnitRange.(iu1, iu2)
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
    get_inputs(model; kwargs...)

Return the input vector corresponding to `model` using the input values in `kwargs`.
"""
function get_inputs(model::AbstractModel; kwargs...)
    y = zeros(number_of_inputs(model))
    set_inputs!(y, model; kwargs...)
    return y
end

"""
    set_inputs!(y, model; kwargs...)

In-place version of [`get_inputs`](@ref)
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

    xs = view.(x, state_indices.(models))

    return ntuple(i->separate_states(xs[i], models[i]), N)
end

"""
    separate_inputs(model, y)

Separate the input vector entries in `y` corresponding to each model in `model`
"""
separate_inputs

separate_inputs(model::NoStateModel, y) = ()

function separate_inputs(models::NTuple{N,AbstractModel}, y) where N

    ys = view.(y, input_indices.(models))

    return ntuple(i->separate_inputs(ys[i], models[i]), N)
end

"""
    separate_parameters(model, p)

Separate the parameter vector entries in `p` corresponding to each model in `model`
and the additional parameters corresponding to the model coupling, if applicable.
"""
separate_parameters

function separate_parameters(models::NTuple{N,AbstractModel}, p) where N

    pmodels = view.(p, parameter_indices.(models))
    pmodels_sep = ntuple(i->separate_parameters(models[i], pmodels[i]), N)

    padd = view(p, additional_parameter_indices(models))
    padd_sep = separate_additional_parameters(models..., padd)

    return (pmodels_sep..., padd_sep)
end

"""
    get_mass_matrix(models)
    get_mass_matrix(models, u, y, p, t)

Calculate the mass matrix for a model or combination of models.
"""
function get_mass_matrix(models::TM, args...; kwargs...) where TM
    return _get_mass_matrix(mass_matrix_type(TM), inplaceness(TM), models,
        args...; kwargs...)
end

# dispatch to an in-place function
function _get_mass_matrix(::Any, ::InPlace, models; kwargs...)
    Nu = number_of_states(models)
    M = zeros(Nu,Nu)
    return get_mass_matrix!(M, models; kwargs...)
end

# dispatch to an in-place function
function _get_mass_matrix(::Any, ::InPlace, models, u, y, p, t; kwargs...)
    Nu = number_of_states(models)
    M = similar(u, Nu, Nu)
    return get_mass_matrix!(M, models, u, y, p, t; kwargs...)
end

# return an empty matrix
function _get_mass_matrix(::Empty, ::OutOfPlace, models, args...; kwargs...)
    return SMatrix{0, 0, Float64}()
end

# return the identity matrix
function _get_mass_matrix(::Zeros, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Nu = number_of_states(TM)

    return zero(SMatrix{Nu, Nu, Float64})
end

# return the identity matrix
function _get_mass_matrix(::Identity, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Nu = number_of_states(TM)

    return SMatrix{Nu, Nu, Float64}(I)
end

# dispatch to the `get_mass_matrix` without arguments
function _get_mass_matrix(::Constant, ::OutOfPlace, models, u, y, p, t; kwargs...)
    return get_mass_matrix(models; kwargs...)
end

# calculate mass matrix for a combination of models
function _get_mass_matrix(::Constant, ::OutOfPlace, models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    M = initialize_mass_matrix(models)

    # calculate input jacobian
    D = input_jacobian(models)

    # calculate input mass matrix
    My = get_coupling_mass_matrix(models)

    return M + D*My
end

# calculate mass matrix for a combination of models
function _get_mass_matrix(::Linear, ::OutOfPlace, models::NTuple{N, AbstractModel},
    u, y, p, t) where N

    # initialize mass matrix
    M = initialize_mass_matrix(models, u, y, p, t)

    # calculate input jacobian
    D = get_input_jacobian(models, u, y, p, t)

    # calculate input mass matrix
    My = get_coupling_mass_matrix(models, u, p, t)

    return M + D*My
end

function initialize_mass_matrix(models::NTuple{N,AbstractModel}) where N
    initialize_static_mass_matrix(models)
end

function initialize_mass_matrix(models::NTuple{N,AbstractModel}, u, y, p, t) where N
    initialize_varying_mass_matrix(models, u, y, p, t)
end

@generated function initialize_static_mass_matrix(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nu = number_of_states.(models.parameters)

    # initialize row names
    Mi = [Symbol("M", i) for i = 1:N]

    # initialize off-diagonal terms
    Mij = [zeros(SMatrix{Nu[i], Nu[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Mi[1]) = vcat(get_mass_matrix(models[1]), $(Mij[1, 2:end]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Mi[i]) = vcat($(Mij[i, 1:i-1]...), get_mass_matrix(models[$i]), $(Mij[i, i+1:end]...))
        end
    end
    expr = quote
        $expr
        M = hcat($(Mi...))
    end

    return expr
end

@generated function initialize_varying_mass_matrix(models::NTuple{N,AbstractModel},
    u, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Mi = [Symbol("M", i) for i = 1:N]

    # initialize off-diagonal terms
    Mij = [zeros(SMatrix{Nu[i], Nu[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Mi[1]) = vcat(get_mass_matrix(models[1], u[$(iu[1])], y[$(iy[1])], p[$(ip[1])], t),
            $(Mij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Mi[i]) = vcat($(Mij[1:i-1, i]...), get_mass_matrix(models[$i], u[$(iu[i])],
                y[$(iy[i])], p[$(ip[i])], t), $(Mij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        M = hcat($(Mi...))
    end

    return expr
end

"""
    get_mass_matrix!(M, models)
    get_mass_matrix!(M, models, u, y, p, t)

In-place version of `get_mass_matrix`.
"""
function get_mass_matrix!(M, models::TM, args...; kwargs...) where TM
    return _get_mass_matrix!(M, mass_matrix_type(TM), inplaceness(TM), models,
        args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_mass_matrix!(M, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return M .= get_mass_matrix(models, args...; kwargs...)
end

# return an empty matrix
function _get_mass_matrix!(M, ::Empty, ::InPlace, models, args...; kwargs...)
    return M
end

# return the identity matrix
function _get_mass_matrix!(M, ::Zeros, ::InPlace, models, args...; kwargs...)
    return M .= 0
end

# return the identity matrix
function _get_mass_matrix!(M, ::Identity, ::InPlace, models, args...; kwargs...)
    M .= 0
    for i = 1:N
        M[i,i] = 1
    end
    return M
end

# dispatch to `get_mass_matrix!` without arguments
function _get_mass_matrix!(M, ::Constant, ::InPlace, models, u, y, p, t; kwargs...)
    return get_mass_matrix!(M, models; kwargs...)
end

# calculate mass matrix for a combination of models
function _get_mass_matrix!(M, ::Constant, ::InPlace, models::NTuple{N, AbstractModel};
    My = similar(M, number_of_inputs(models), number_of_states(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)

    # calculate input mass matrix
    get_coupling_mass_matrix!(My, models)

    # calculate mass matrix
    for i = 1:length(models)
        D = get_input_jacobian(models[i])
        for j = 1:length(models)
            Mij = view(M, iu[i], iu[j])
            Myij = view(My, iy[i], iu[j])
            if i == j
                get_mass_matrix!(Mij, models)
                mul!(Mij, D, Myij, 1, 1)
            else
                mul!(Mij, D, Myij)
            end
        end
    end

    return M
end

# calculate mass matrix for a combination of models
function _get_mass_matrix!(M, ::Linear, ::InPlace,
    models::NTuple{N, AbstractModel}, u, y, p, t; My = similar(M, number_of_inputs(models),
    number_of_states(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input mass matrix
    get_coupling_mass_matrix!(My, models, u, p, t)

    # calculate mass matrix
    for i = 1:N
        ui = view(u, iu[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        D = get_input_jacobian(models[i], ui, yi, pi, t)
        for j = 1:N
            Mij = view(M, iu[i], iu[j])
            Myij = view(My, iy[i], iu[j])
            if i == j
                get_mass_matrix!(Mij, models[i], ui, yi, pi, t)
                mul!(Mij, D, Myij, 1, 1)
            else
                mul!(Mij, D, Myij)
            end
        end
    end

    return M
end

"""
    get_lhs(model, u, y, p, t)

Return the left hand side of the governing differential equations for the model.
This method is used for testing the mass matrix associated with each model.
"""
get_lhs(models, u, y, p, t)

"""
    get_rates(models, u, y, p, t)

Calculate the (mass matrix multiplied) state rates for the specified model or
models.
"""
function get_rates(models::T, u, y, p, t) where T
    return _get_rates(inplaceness(T), models, u, y, p, t)
end

function get_rates(models::T, args...) where T<:NoStateModel
    return SVector{0,Float64}()
end

# dispatch to an in-place function
function _get_rates(::InPlace, models, u, y, p, t)
    TF = promote_type(eltype(u), eltype(y), eltype(p), typeof(t))
    du = similar(u, TF)
    get_rates!(du, models, u, y, p, t)
    return du
end

# calculate state rates for a combination of models
function _get_rates(::OutOfPlace, models::NTuple{N,AbstractModel}, u, y, p, t) where N
    return get_model_rates(models, u, y, p, t)
end

@generated function get_model_rates(models::NTuple{N,AbstractModel}, u, y, p, t) where N
    # get indices of state, input, and parameter vectors for each model
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize state variable names
    dui = [Symbol("du", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(dui[i]) = get_rates(models[$i], u[$(iu[i])], y[$(iy[i])], p[$(ip[i])], t)
        end
    end

    expr = quote
        $expr
        du = vcat($(dui...))
    end

    return expr
end

"""
    get_rates!(du, models, u, y, p, t)

In-place version of [`get_rates`](@ref)
"""
function get_rates!(du, models::T, u, y, p, t) where T
    return _get_rates!(du, inplaceness(T), models, u, y, p, t)
end

# dispatch to an out-of-place function
function _get_rates!(du, ::OutOfPlace, models, u, y, p, t)
    return du .= get_rates(models, u, y, p, t)
end

# calculate state rates for a combination of models
function _get_rates!(du, ::InPlace, models::NTuple{N,AbstractModel}, u, y, p, t) where N

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    for i = 1:length(models)
        vdu = view(du, iu[i])
        vu = view(u, iu[i])
        vy = view(y, iy[i])
        vp = view(p, ip[i])
        get_rates!(vdu, models[i], vu, vy, vp, t)
    end

    return du
end

"""
    get_state_jacobian(models, u, y, p, t)

Calculate the jacobian with respect to the state variables for the specified
models.
"""
function get_state_jacobian(models::T, u, y, p, t) where T
    return _get_state_jacobian(state_jacobian_type(T), inplaceness(T), models,
        u, y, p, t)
end

# dispatch to an in-place function
function _get_state_jacobian(::Any, ::InPlace, models)
    Nu = number_of_states(models)
    M = zeros(Nu,Nu)
    return get_state_jacobian!(M, models)
end

# dispatch to an in-place function
function _get_state_jacobian(::Any, ::InPlace, models, u, y, p, t)
    TF = promote_type(eltype(u), eltype(y), eltype(p), typeof(t))
    Nu = number_of_states(models)
    M = zeros(TF, Nu, Nu)
    return get_state_jacobian!(M, models, u, y, p, t)
end

# return an empty matrix
function _get_state_jacobian(::Empty, ::OutOfPlace, models, args...)
    return SMatrix{0, 0, Float64}()
end

# return a zero matrix
function _get_state_jacobian(::Zeros, ::OutOfPlace, models::TM, args...) where TM
    Nu = number_of_states(TM)
    return zero(SMatrix{Nu, Nu, Float64})
end

# return the identity matrix
function _get_state_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nu = number_of_states(TM)
    return SMatrix{Nu, Nu, Float64}(I)
end

# dispatch to the `get_state_jacobian` without arguments
function _get_state_jacobian(::Constant, ::OutOfPlace, models, u, y, p, t)
    return get_state_jacobian(models)
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace, models,
    u, y, p, t)

    f = u -> get_rates(models, u, y, p, t)

    return ForwardDiff.jacobian(f, u)
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    J = initialize_state_jacobian(models)

    # calculate input jacobian
    D = input_jacobian(models)

    # calculate input jacobian
    Jy = get_coupling_state_jacobian(models)

    return J + D*Jy
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, u, y, p, t) where N

    # initialize mass matrix
    J = initialize_state_jacobian(models, u, y, p, t)

    # calculate input jacobian
    D = get_input_jacobian(models, u, y, p, t)

    # calculate input jacobian
    Jy = get_coupling_state_jacobian(models, u, p, t)

    return J + D*Jy
end

function initialize_state_jacobian(models::NTuple{N,AbstractModel}) where N
    initialize_static_state_jacobian(models)
end

function initialize_state_jacobian(models::NTuple{N,AbstractModel}, u, y, p, t) where N
    initialize_varying_state_jacobian(models, u, y, p, t)
end

@generated function initialize_static_state_jacobian(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nu = number_of_states.(models.parameters)

    # initialize row names
    Ji = [Symbol("J", i) for i = 1:N]

    # initialize off-diagonal terms
    Jij = [zeros(SMatrix{Nu[i], Nu[j], Float64}) for i = 1:N, j = 1:N]

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

@generated function initialize_varying_state_jacobian(models::NTuple{N,AbstractModel},
    u, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)

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
    Jij = [zeros(SMatrix{Nu[i], Nu[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Ji[1]) = vcat(get_state_jacobian(models[1], u[$(iu[1])], y[$(iy[1])],
            p[$(ip[1])], t), $(Jij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Ji[i]) = vcat($(Jij[1:i-1, i]...), get_state_jacobian(models[$i],
                u[$(iu[i])], y[$(iy[i])], p[$(ip[i])], t), $(Jij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        J = hcat($(Ji...))
    end

    return expr
end

"""
    get_state_jacobian!(J, models, u, y, p, t)

In-place version of [`get_state_jacobian`](@ref)
"""
function get_state_jacobian!(J, models::T, args...; kwargs...) where T
    return _get_state_jacobian!(J, state_jacobian_type(T), inplaceness(T),
        models, args...; kwargs...)
end

# dispatch to an out-of-place function
function _get_state_jacobian!(J, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return J .= get_state_jacobian(models, args...; kwargs...)
end

# return an empty matrix
function _get_state_jacobian!(J, ::Empty, ::InPlace, models, args...; kwargs...)
    return J
end

# return the identity matrix
function _get_state_jacobian!(J, ::Zeros, ::InPlace, models, args...; kwargs...)
    return J .= 0
end

# return the identity matrix
function _get_state_jacobian!(J, ::Identity, ::InPlace, models, args...; kwargs...)
    J .= 0
    for i = 1:N
        J[i,i] = 1
    end
    return J
end

# dispatch to `get_state_jacobian!` without arguments
function _get_state_jacobian!(J, ::Constant, ::OutOfPlace, models, u, y, p, t; kwargs...)
    return get_state_jacobian!(J, models; kwargs...)
end

# calculate state jacobian for a combination of models
function _get_state_jacobian!(J, ::Constant, ::InPlace, models::NTuple{N,AbstractModel};
    Jy = zeros(number_of_inputs(models), number_of_states(models))) where N

    # get dimensions
    Nu = number_of_states.(models)
    Ny = number_of_inputs.(models)
    Np = number_of_parameters.(models)

    # state variable indices
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = UnitRange.(iu1, iu2)

    # input variable indices
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = UnitRange.(iy1, iy2)

    # get input jacobian
    get_coupling_state_jacobian!(Jy, models)

    # calculate jacobian
    for i = 1:length(models)
        D = get_input_jacobian(models[i])
        for j = 1:length(models)
            Jij = view(J, iu[i], iu[j])
            Jyij = view(Jy, iy[i], iu[j])
            if i == j
                get_state_jacobian!(Jij, models[i])
                mul!(Jij, D, Jyij, 1, 1)
            else
                mul!(Jij, D, Jyij)
            end
        end
    end

    return J
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian!(J, ::Union{Linear, Nonlinear}, ::InPlace, model,
    u, y, p, t)

    f = u -> get_rates(models, u, y, p, t)

    return copyto!(J, ForwardDiff.jacobian(f, u))
end

# calculate state jacobian for a combination of models
function _get_state_jacobian!(J, ::Union{Linear, Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, u, y, p, t;
    Jy = zeros(number_of_inputs(models), number_of_states(models))) where N

    # get dimensions
    Nu = number_of_states.(models)
    Ny = number_of_inputs.(models)
    Np = number_of_parameters.(models)

    # state variable indices
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = UnitRange.(iu1, iu2)

    # input variable indices
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = UnitRange.(iy1, iy2)

    # parameter indices
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = UnitRange.(ip1, ip2)

    # get input jacobian
    get_coupling_state_jacobian!(Jy, models, u, p, t)

    # calculate jacobian
    for i = 1:length(models)
        ui = view(u, iu[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        D = get_input_jacobian(models[i], ui, yi, pi, t)
        for j = 1:length(models)
            Jij = view(J, iu[i], iu[j])
            Jyij = view(Jy, iy[i], iu[j])
            if i == j
                get_state_jacobian!(Jij, models[i], ui, yi, pi, t)
                mul!(Jij, D, Jyij, 1, 1)
            else
                mul!(Jij, D, Jyij)
            end
        end
    end

    return J
end

"""
    get_input_jacobian(models)
    get_input_jacobian(models, u, y, p, t)

Calculate the jacobian with respect to the inputs for the specified model or models.
"""
function get_input_jacobian(models::TM, args...) where TM
    return _get_input_jacobian(input_jacobian_type(TM), models, args...)
end

# return an empty matrix
function _get_input_jacobian(::Empty, models::TM, args...) where TM
    Ny = number_of_inputs(TM)
    return SMatrix{Ny, 0, Float64}()
end

# return a zero-valued matrix
function _get_input_jacobian(::Zeros, models, args...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    return FillMap(0, Ny, Nu)
end

# return an identity matrix
function _get_input_jacobian(::Identity, models, args...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    @assert Ny == Nu
    return LinearMap(I, Ny)
end

# dispatch to the `get_input_jacobian` function without arguments
function _get_input_jacobian(::Constant, models, u, y, p, t)
    return get_input_jacobian(models)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Constant, models::NTuple{N,AbstractModel}) where N
    static_input_jacobian(models)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian(::Union{Linear, Nonlinear}, model, u, y, p, t)

    f = y -> get_rates(model, u, y, p, t)

    return ForwardDiff.jacobian(f, y)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Union{Linear, Nonlinear},
    models::NTuple{N,AbstractModel}, u, y, p, t) where N

    varying_input_jacobian(models, u, y, p, t)
end

@generated function static_input_jacobian(models::NTuple{N,AbstractModel}) where N

    # get number of state variables in each model
    Nu = number_of_states.(models.parameters)
    Ny = number_of_inputs.(models.parameters)

    # initialize row names
    Di = [Symbol("D", i) for i = 1:N]

    # initialize off-diagonal terms
    Dij = [zeros(SMatrix{Nu[i], Ny[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Di[1]) = vcat(get_input_jacobian(models[1]), $(Dij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Di[i]) = vcat($(Dij[1:i-1, i]...), get_input_jacobian(models[$i]),
                $(Dij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        D = hcat($(Di...))
    end

    return expr
end

@generated function varying_input_jacobian(models::NTuple{N,AbstractModel},
    u, y, p, t) where N

    # get indices of state, input, and parameter vectors for each model
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)

    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)

    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)

    # initialize row names
    Di = [Symbol("D", i) for i = 1:N]

    # initialize off-diagonal terms
    Dij = [zeros(SMatrix{Nu[i], Ny[j], Float64}) for i = 1:N, j = 1:N]

    # construct all columns
    expr = quote
        $(Di[1]) = vcat(get_input_jacobian(models[1], u[$(iu[1])], y[$(iy[1])], p[$(ip[1])], t),
            $(Dij[2:end, 1]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Di[i]) = vcat($(Dij[1:i-1, i]...), get_input_jacobian(models[$i], u[$(iu[i])],
                y[$(iy[i])], p[$(ip[i])], t), $(Dij[i+1:end, i]...))
        end
    end
    expr = quote
        $expr
        D = hcat($(Di...))
    end

    return expr
end

"""
    get_coupling_mass_matrix(models)
    get_coupling_mass_matrix(models, u, p, t)

Calculate the input function mass matrix for the specified combination of models.
"""
get_coupling_mass_matrix

function get_coupling_mass_matrix(models::TM; kwargs...) where TM
    return _get_coupling_mass_matrix(coupling_mass_matrix_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models; kwargs...)
end

function get_coupling_mass_matrix(models::TM, u, p, t; kwargs...) where TM
    return _get_coupling_mass_matrix(coupling_mass_matrix_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an in-place function
function _get_coupling_mass_matrix(::Any, ::InPlace, models; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = zeros(Ny,Nu)
    return get_coupling_mass_matrix!(M, models; kwargs...)
end

# dispatch to an in-place function
function _get_coupling_mass_matrix(::Any, ::InPlace, models, u, p, t; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = similar(u, Ny, Nu)
    return get_coupling_mass_matrix!(M, models, u, p, t; kwargs...)
end

# return an empty matrix
function _get_coupling_mass_matrix(::Empty, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)

    return SMatrix{Ny, 0, Float64}()
end

# return a matrix of zeros
function _get_coupling_mass_matrix(::Zeros, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)
    Nu = number_of_states(TM)

    return zeros(SMatrix{Ny, Nu, Float64})
end

# dispatch to the `get_coupling_mass_matrix` function without arguments
function _get_coupling_mass_matrix(::Constant, ::OutOfPlace, models, u, p, t; kwargs...)
    return get_coupling_mass_matrix(models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_mass_matrix(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    return get_coupling_mass_matrix(models...; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_mass_matrix(::Linear, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    return get_coupling_mass_matrix(models..., u, p, t; kwargs...)
end

"""
    get_coupling_mass_matrix!(My, models)
    get_coupling_mass_matrix!(My, models, u, p, t)

In-place version of [`get_coupling_mass_matrix`](@ref).
"""
get_coupling_mass_matrix!

function get_coupling_mass_matrix!(My, models::TM; kwargs...) where TM
    return _get_coupling_mass_matrix!(My, coupling_mass_matrix_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models; kwargs...)
end

function get_coupling_mass_matrix!(My, models::TM, u, p, t; kwargs...) where TM
    return _get_coupling_mass_matrix!(My, coupling_mass_matrix_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_coupling_mass_matrix!(My, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return My .= get_coupling_mass_matrix(models, args...; kwargs...)
end

# return an empty matrix
function _get_coupling_mass_matrix!(My, ::Empty, ::InPlace, models, args...; kwargs...)
    return My
end

# return a matrix of zeros
function _get_coupling_mass_matrix!(My, ::Zeros, ::InPlace, models, args...; kwargs...)
    return My .= 0
end

# dispatch to `get_coupling_mass_matrix!` without arguments
function _get_coupling_mass_matrix!(My, ::Constant, ::InPlace, models, u, p, t; kwargs...)
    return get_coupling_mass_matrix!(My, models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_mass_matrix!(My, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    get_coupling_mass_matrix!(My, models...; kwargs...)

    return My
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_mass_matrix!(My, ::Linear, ::InPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    get_coupling_mass_matrix!(My, models..., u, p, t; kwargs...)

    return My
end

"""
    get_inputs_using_state_rates(models, u, p, t)

Calculate the portion of the inputs which are dependent on the state rates.  This
function is used to test the mass matrices associated with each input mass matrix.
"""
function get_inputs_using_state_rates(models::NTuple{N,AbstractModel}, du, u, p, t) where N
    return get_inputs_using_state_rates(models..., du, u, p, t)
end

"""
    get_inputs(models, u, p, t)

Calculate the inputs to the specified combination of models.
"""
function get_inputs(models::TM, u, p, t) where TM
    return _get_inputs(coupling_inplaceness(TM.parameters...), models, u, p, t)
end

# dispatch to an in-place function
function _get_inputs(::InPlace, models, u, p, t)
    Ny = number_of_inputs(models)
    y = similar(u, Ny)
    get_inputs!(y, models, u, p, t)
    return y
end

# dispatch to the user-provided function for the specific combination of models
function _get_inputs(::OutOfPlace, models::NTuple{N,AbstractModel}, u, p, t) where N
    return get_inputs(models..., u, p, t)
end

"""
    get_inputs!(y, models::NTuple{N,AbstractModel}, u, p, t) where N

In-place version of [`get_inputs`](@ref)
"""
function get_inputs!(y, models::T, u, p, t) where T
    return _get_inputs!(y, coupling_inplaceness(T.parameters...), models, u, p, t)
end

# dispatch to an out-of-place function
function _get_inputs!(y, ::OutOfPlace, models, u, p, t)
    return y .= get_inputs(models, u, p, t)
end

# dispatch to the user-provided function for the specific combination of models
function _get_inputs!(y, ::InPlace, models::NTuple{N,AbstractModel}, u, p, t) where N
    return get_inputs!(y, models..., u, p, t)
end

"""
    get_coupling_state_jacobian(models::NTuple{N,AbstractModel}, u, p, t) where N

Calculate the jacobian of the input function with respect to the state variables
for the specified combination of models.
"""
get_coupling_state_jacobian

function get_coupling_state_jacobian(models::TM; kwargs...) where TM
    return _get_coupling_state_jacobian(coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models; kwargs...)
end

function get_coupling_state_jacobian(models::TM, u, p, t; kwargs...) where TM
    return _get_coupling_state_jacobian(coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_state_jacobian(args...)

    models = args[1:end-3]
    u = args[end-2]
    p = args[end-1]
    t = args[end]

    f = u -> get_inputs(models, u, p, t)

    return ForwardDiff.jacobian(f, u)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::Any, ::InPlace, models; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = zeros(Ny,Nu)
    return get_coupling_state_jacobian!(M, models; kwargs...)
end

# dispatch to an in-place function
function _get_coupling_state_jacobian(::Any, ::InPlace, models, u, p, t; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = similar(u, Ny, Nu)
    return get_coupling_state_jacobian!(M, models, u, p, t; kwargs...)
end

# return an empty matrix
function _get_coupling_state_jacobian(::Empty, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)

    return SMatrix{Ny, 0, Float64}()
end

# return a matrix of zeros
function _get_coupling_state_jacobian(::Zeros, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)
    Nu = number_of_states(TM)

    return zeros(SMatrix{Ny, Nu, Float64})
end

# dispatch to the `get_coupling_state_jacobian` function without arguments
function _get_coupling_state_jacobian(::Constant, ::OutOfPlace, models, u, p, t; kwargs...)
    return get_coupling_state_jacobian(models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    return get_coupling_state_jacobian(models...; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    return get_coupling_state_jacobian(models..., u, p, t; kwargs...)
end

"""
    get_coupling_state_jacobian!(J, models, u, p, t)

In-place version of [`get_coupling_state_jacobian`](@ref)
"""
get_coupling_state_jacobian!

function get_coupling_state_jacobian!(Jy, models::TM; kwargs...) where TM
    return _get_coupling_state_jacobian!(Jy, coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models; kwargs...)
end

function get_coupling_state_jacobian!(Jy, models::TM, u, p, t; kwargs...) where TM
    return _get_coupling_state_jacobian!(Jy, coupling_state_jacobian_type(TM.parameters...),
        coupling_inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# use automatic differentiation if jacobian is not defined
function get_coupling_state_jacobian!(Jy, args...)

    models = args[1:end-3]
    u = args[end-2]
    p = args[end-1]
    t = args[end]

    f = u -> get_inputs(models, u, p, t)

    return copyto!(Jy, ForwardDiff.jacobian(f, u))
end

# dispatch to an out-of-place function
function _get_coupling_state_jacobian!(Jy, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return Jy .= get_coupling_state_jacobian(models, args...; kwargs...)
end

# return an empty matrix
function _get_coupling_state_jacobian!(Jy, ::Empty, ::InPlace, models, args...; kwargs...)
    return Jy
end

# return a matrix of zeros
function _get_coupling_state_jacobian!(Jy, ::Zeros, ::InPlace, models, args...; kwargs...)
    return Jy .= 0
end

# dispatch to `get_coupling_state_jacobian!` without arguments
function _get_coupling_state_jacobian!(Jy, ::Constant, ::InPlace, models, u, p, t; kwargs...)
    return get_coupling_state_jacobian!(Jy, models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian!(Jy, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    get_coupling_state_jacobian!(Jy, models...; kwargs...)

    return Jy
end

# dispatch to the user-provided function for the specific combination of models
function _get_coupling_state_jacobian!(Jy, ::Union{Linear, Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    get_coupling_state_jacobian!(Jy, models..., u, p, t; kwargs...)

    return Jy
end

"""
    get_eigen(model::TM, x, y, p, t; kwargs...)

Return the eigenvalues, left eigenvector matrix, and right eigenvector matrix
corresponding to the model for state variables `x`, inputs `y`, parameters `p`,
and time `t`.

For in-place models, the number of eigenvalues to compute may be specified using
the `nev` keyword argument.
"""
function get_eigen(model::TM, x, y, p, t; kwargs...) where TM
    return _get_eigen(inplaceness(TM), model, x, y, p, t; kwargs...)
end

function _get_eigen(::OutOfPlace, model, x, y, p, t)
    M = Array(get_mass_matrix(model, x, y, p, t)) # mass matrix
    K = Array(get_state_jacobian(model, x, y, p, t)) # jacobian
    E = eigen(K, M) # eigenvalue decomposition
    λ = eigvals(E) # eigenvalues
    V = eigvecs(E) # right eigenvector matrix
    U = I/(M*V) # left eigenvector matrix
    return λ, U, V
end

function _get_eigen(::InPlace, model, x, y, p, t; nev=20)

    # calculate the mass matrix corresponding to steady state operating conditions
    M = get_mass_matrix(model, x, y, p, t)

    # calculate the jacobian corresponding to steady state operating conditions
    K = get_state_jacobian(model, x, y, p, t)

    # construct linear map
    T = promote_type(eltype(K), eltype(M))
    nx = size(K, 1)
    Kfact = lu(K)
    f! = (b, x) -> ldiv!(b, Kfact, M*x)
    fc! = (b, x) -> mul!(b, M', Kfact'\x)
    A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V = partialeigen(partialschur(A; nev=min(nx,nev), which=LM())[1])

    # sort eigenvalues by magnitude
    perm = sortperm(λ, by=(λ)->(abs(λ),imag(λ)), rev=true)
    λ .= λ[perm]
    V .= V[:,perm]

    # eigenvalues are actually 1/λ, no modification necessary for eigenvectors
    λ .= 1 ./ λ

    # also return left eigenvectors
    U = GXBeam.left_eigenvectors(K, -M, λ, V)

    return λ, U, V
end

"""
    get_ode(model)

Construct an ODEFunction corresponding to the specified model or models which
may be solved using DifferentialEquations.
"""
function get_ode(model::TM) where TM
    return _get_ode(mass_matrix_type(TM), inplaceness(TM), model)
end

function _get_ode(::Identity, ::OutOfPlace, model::AbstractModel)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (u, p, t) -> get_rates(models, u, p[iy], p[ip], t)

    jac = (u, p, t) -> get_state_jacobian(models, u, p[iy], p[ip], t)

    return ODEFunction{false}(f; jac)
end

function _get_ode(::Constant, ::OutOfPlace, model::AbstractModel)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (u, p, t) -> get_rates(models, u, p[iy], p[ip], t)

    mass_matrix = get_mass_matrix(models)

    jac = (u, p, t) -> get_state_jacobian(models, u, p[iy], p[ip], t)

    return ODEFunction{false}(f; mass_matrix, jac)
end

function _get_ode(::Linear, ::OutOfPlace, model::AbstractModel)

    Np = number_of_parameters(model)
    Ny = number_of_inputs(model)

    ip = SVector{Np}(1:Np)
    iy = SVector{Ny}(Np+1:Np+Ny)

    f = (du, u, p, t) -> get_rates!(du, models, u, p[iy], p[ip], t)

    M = zeros(Nu, Nu)
    update_func = (M, u, p, t) -> get_mass_matrix!(M, models, p[iy], p[ip], t)
    mass_matrix = DiffEqArrayOperator(M; update_func)

    jac = (J, u, p, t) -> get_state_jacobian!(J, models, u, p[iy], p[ip], t)

    return ODEFunction{true}(f; mass_matrix, jac)
end

function _get_ode(::Identity, ::InPlace, model::AbstractModel)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    # construct state rate function
    f = (du, u, p, t) ->  get_rates!(du, model, u, view(p, iy), view(p, ip), t)

    # construct jacobian function
    jac = (J, u, p, t) -> get_state_jacobian!(J, model, u, view(p, iy), view(p, ip), t)

    # construct and return an ODEFunction
    return ODEFunction{true}(f; jac)
end

function _get_ode(::Constant, ::InPlace, model::AbstractModel)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    # construct state rate function
    f = (du, u, p, t) -> get_rates!(du, model, u, view(p, iy), view(p, ip), t)

    # construct mass matrix
    mass_matrix = get_mass_matrix(model)

    # construct jacobian function
    jac = (J, u, p, t) -> get_state_jacobian!(J, model, u, view(p, iy), view(p, ip), t)

    # construct and return an ODEFunction
    return ODEFunction{true}(f; mass_matrix, jac)
end

function _get_ode(::Linear, ::InPlace, model::AbstractModel)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    ip = 1:Np
    iy = Np+1:Np+Ny

    # construct state rate function
    f = (du, u, p, t) -> get_rates!(du, model, u, view(p, iy), view(p, ip), t)

    # construct mass matrix
    M = zeros(Nu, Nu)
    update_func = (M, u, p, t) -> get_mass_matrix!(M, model, u, view(p, iy), view(p, ip), t)
    mass_matrix = DiffEqArrayOperator(M; update_func)

    # construct jacobian function
    jac = (J, u, p, t) -> get_state_jacobian!(J, model, u, view(p, iy), view(p, ip), t)

    # construct and return an ODEFunction
    return ODEFunction{true}(f; mass_matrix, jac)
end

function _get_ode(::Identity, ::OutOfPlace, model::Tuple)

    fy = (u, p, t) -> get_inputs(models, u, p, t)

    f = (u, p, t) -> get_rates(models, u, fy(u, p, t), p, t)

    jac = (u, p, t) -> get_state_jacobian(models, u, fy(u, p, t), p, t)

    return ODEFunction{false}(f; jac)
end

function _get_ode(::Constant, ::OutOfPlace, model::Tuple)

    fy = (u, p, t) -> get_inputs(models, u, p, t)

    f = (u, p, t) -> get_rates(models, u, fy(u, p, t), p, t)

    mass_matrix = get_mass_matrix(models)

    jac = (u, p, t) -> get_state_jacobian(models, u, fy(u, p, t), p, t)

    return ODEFunction{false}(f; mass_matrix, jac)
end

function _get_ode(::Linear, ::OutOfPlace, model::Tuple)

    Nu = number_of_states(model)

    fy = (u, p, t) -> get_inputs(model, u, p, t)

    f = (du, u, p, t) -> get_rates!(du, model, u, fy(u, p, t), p, t)

    M = zeros(Nu, Nu)
    update_func = (M, u, p, t) -> get_mass_matrix!(M, model, u, fy(u, p, t), p, t)
    mass_matrix = DiffEqArrayOperator(M; update_func)

    jac = (J, u, p, t) -> get_state_jacobian!(J, model, u, fy(u, p, t), p, t)

    return ODEFunction{true}(f; mass_matrix, jac)
end

function _get_ode(::Identity, ::InPlace, model::Tuple)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # construct state rate function
    ucache = fill(NaN, Nu)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)
    f = (du, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate state rates
        get_rates!(du, model, u, y, p, t)
    end

    # construct jacobian function
    jac = (J, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate jacobian
        get_state_jacobian!(J, model, u, y, p, t)
    end

    # construct and return an ODEFunction
    return ODEFunction{true}(f; jac)
end

function _get_ode(::Constant, ::InPlace, model::Tuple)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # construct state rate function
    ucache = fill(NaN, Nu)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)
    f = (du, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate state rates
        get_rates!(du, model, u, y, p, t)
    end

    # construct mass matrix
    mass_matrix = get_mass_matrix(model)

    # construct jacobian function
    jac = (J, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate jacobian
        get_state_jacobian!(J, model, u, y, p, t)
    end

    # construct and return an ODEFunction
    return ODEFunction{true}(f; mass_matrix, jac)
end

function _get_ode(::Linear, ::InPlace, model::Tuple)

    # problem dimensions
    Nu = number_of_states(model)
    Ny = number_of_inputs(model)
    Np = number_of_parameters(model)

    # construct state rate function
    ucache = fill(NaN, Nu)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)
    f = (du, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate state rates
        get_rates!(du, model, u, y, p, t)
    end

    # construct mass matrix
    M = zeros(Nu, Nu)
    My_cache = zeros(Ny, Nu)
    update_func = (M, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(M) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                # store current input arguments
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # update type of `My`
        My = convert(typeof(M), My_cache)
        # calculate inputs
        get_mass_matrix!(M, model, u, y, p, t; My)
    end
    mass_matrix = DiffEqArrayOperator(M; update_func)

    # construct jacobian function
    jac = (J, u, p, t) -> begin
        # check if we can use the cache variables (no custom types)
        if eltype(du) <: Float64
            # update the cache variables
            if (u != ucache) && (p != pcache) && (t != tcache[])
                ucache .= u
                get_inputs!(ycache, model, u, p, t)
                pcache .= p
                tcache .= t
            end
            # use the cached model inputs
            y = ycache
        else
            # calculate model inputs (out-of-place to accomodate custom type)
            y = get_inputs(model, u, p, t)
        end
        # calculate jacobian
        get_state_jacobian!(J, model, u, y, p, t)
    end

    # construct and return an ODEFunction
    return ODEFunction{true}(f; mass_matrix, jac)
end

# plotting recipes
@recipe function f(models::AbstractModel, x, p, t) where N

    Np = number_of_parameters(models)
    Ny = number_of_inputs(models)

    ip = 1:Np
    iy = Np+1:Np+Ny

    return models..., x, p[iy], p[ip], t
end

@recipe function f(models::NTuple{N,AbstractModel}, x, p, t) where N

    y = get_inputs(models, x, p, t)

    return models..., x, y, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, sol) where N

    it = sol.tslocation

    x = sol.u[it]
    p = sol.prob.p
    t = sol.t[it]

    y = get_inputs(models, x, p, t)

    return models..., x, y, p, t
end

@recipe function f(models::NTuple{N,AbstractModel}, sol, t) where N

    x = sol(t)
    p = sol.prob.p

    y = get_inputs(models, x, p, t)

    return models..., x, y, p, t
end
