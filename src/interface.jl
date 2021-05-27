"""
    AbstractModel

Supertype for all aerodynamic, structural, and dynamics models.
"""
abstract type AbstractModel end

"""
    AerodynamicModel

Supertype for all aerodynamic models
"""
abstract type AerodynamicModel <: AbstractModel end

"""
    StructuralModel

Supertype for all structural models
"""
abstract type StructuralModel <: AbstractModel end

"""
    DynamicsModel

Supertype for all dynamics models
"""
abstract type DynamicsModel <: AbstractModel end

"""
    number_of_states(models)

Return the total number of states corresponding to the model or models.
"""
number_of_states

number_of_states(model::AbstractModel) = number_of_states(typeof(model))

function number_of_states(models::NTuple{N,AbstractModel}) where N
    sum(number_of_states.(models))
end

function number_of_states(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_states.((T.parameters...)))
end

"""
    number_of_inputs(models)

Return the total number of inputs corresponding to the model or models.
"""
number_of_inputs

number_of_inputs(model::AbstractModel) = number_of_inputs(typeof(model))

function number_of_inputs(models::NTuple{N,AbstractModel}) where N
    sum(number_of_inputs.(models))
end

function number_of_inputs(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_inputs.((T.parameters...)))
end

"""
    number_of_parameters(models)

Return the total number of parameters corresponding to the model or models.
"""
number_of_parameters

number_of_parameters(model::AbstractModel) = number_of_parameters(typeof(model))

function number_of_parameters(models::NTuple{N,AbstractModel}) where N
    sum(number_of_parameters.(models))
end

function number_of_parameters(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_parameters.((T.parameters...)))
end

function state_indices(models::NTuple{N,AbstractModel}) where N
    # input indices
    Nu = number_of_states.(typeof.(models))
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    return UnitRange.(iu1, iu2)
end

@generated function static_state_indices(models::NTuple{N,AbstractModel}) where N
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)
    return :($iu)
end

function input_indices(models::NTuple{N,AbstractModel}) where N
    # input indices
    Ny = number_of_inputs.(typeof.(models))
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    return UnitRange.(iy1, iy2)
end

@generated function static_input_indices(models::NTuple{N,AbstractModel}) where N
    Ny = number_of_inputs.(models.parameters)
    iy2 = cumsum(Ny)
    iy1 = iy2 .- Ny .+ 1
    iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), N)
    return :($iy)
end

function parameter_indices(models::NTuple{N,AbstractModel}) where N
    # parameter indices
    Np = number_of_parameters.(typeof.(models))
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    return UnitRange.(ip1, ip2)
end

@generated function static_parameter_indices(models::NTuple{N,AbstractModel}) where N
    Np = number_of_parameters.(models.parameters)
    ip2 = cumsum(Np)
    ip1 = ip2 .- Np .+ 1
    ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), N)
    return :($ip)
end

"""
    isinplace(models)

Return a flag indicating whether the state rates corresponding to the model or
combination of models are calculated in-place.

This property should be inferrable from the model type.
"""
isinplace

isinplace(model::AbstractModel) = isinplace(typeof(model))

function isinplace(models::NTuple{N,AbstractModel}) where N
    inplace_inputs(models...) || any(isinplace.(models))
end

function isinplace(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    inplace_inputs(model_types...) || any(isinplace.(model_types))
end

"""
    has_mass_matrix(models)

Return a flag indicating whether the governing differential equations for the
model or combination or models has a mass matrix.

This property should be inferrable from the model type.
"""
has_mass_matrix

has_mass_matrix(model::AbstractModel) = has_mass_matrix(typeof(model))

function has_mass_matrix(models::NTuple{N,AbstractModel}) where N
    has_input_mass_matrix(models) || any(has_mass_matrix.(models))
end

function has_mass_matrix(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    has_input_mass_matrix(model_types...) || any(has_mass_matrix.(model_types))
end

"""
    constant_mass_matrix(models)

Return a flag indicating whether the mass matrix for the model or combination of
models is constant.

This property should be inferrable from the model type.
"""
constant_mass_matrix

function constant_mass_matrix(model::AbstractModel)
    constant_matrix_state(mass_matrix_state(typeof(model)))
end

function constant_mass_matrix(models::NTuple{N,AbstractModel}) where N
    return constant_input_mass_matrix(models...) &&
        all(constant_input_jacobian.(models)) &&
        all(constant_mass_matrix.(models))
end

function constant_mass_matrix(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    return constant_input_mass_matrix(model_types...) &&
        all(constant_input_jacobian.(model_types)) &&
        all(constant_mass_matrix.(model_types))
end

"""
    linear_input_dependence(models)

Return a flag indicating whether the state rates corresponding to the model or
combination of models has a linear dependence on the inputs
"""
linear_input_dependence

linear_input_dependence(model::AbstractModel) = linear_input_dependence(typeof(model))

function linear_input_dependence(models::NTuple{N,AbstractModel}) where N
    return all(linear_input_dependence.(models))
end

function linear_input_dependence(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    return all(linear_input_dependence.(model_types))
end

"""
    defined_state_jacobian(models)

Return a flag indicating whether the jacobian of the state rates with respect
to the states is manually defined.
"""
defined_state_jacobian

defined_state_jacobian(model::AbstractModel) = defined_state_jacobian(typeof(model))

function defined_state_jacobian(models::NTuple{N,AbstractModel}) where N
    return defined_input_state_jacobian(models...) &&
        all(defined_state_jacobian.(models)) &&
        all(defined_input_jacobian.(models))
end

function defined_state_jacobian(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    model_types = (T.parameters...,)
    return defined_input_state_jacobian(model_types...) &&
        all(defined_state_jacobian.(model_types)) &&
        all(defined_input_jacobian.(model_types))
end

"""
    defined_input_jacobian(model)

Return a flag indicating whether the jacobian of the state rates with respect
to the inputs is defined. Defaults to false.
"""
defined_input_jacobian

defined_input_jacobian(model::AbstractModel) = defined_input_jacobian(typeof(model))

defined_input_jacobian(::Type{T}) where T<:AbstractModel = false

"""
    constant_input_jacobian(model)

Return a flag indicating whether the jacobian of the state rates with respect
to the inputs is constant.  Defaults to false.
"""
function constant_input_jacobian(model::AbstractModel)
    constant_matrix_state(input_jacobian_state(model))
end

constant_input_jacobian(::Type{T}) where T<:AbstractModel = false

"""
    get_mass_matrix(models)
    get_mass_matrix(models, u, y, p, t)

Calculate the mass matrix for a model or combination of models.
"""
get_mass_matrix

function get_mass_matrix(model::AbstractModel, u, y, p, t)
    return get_mass_matrix(mass_matrix_state(model), inplaceness(model), model,
        u, y, p, t)
end

function get_mass_matrix(models::NTuple{N,AbstractModel}, args...) where N
    return get_mass_matrix(mass_matrix_state(models), inplaceness(models), models, args...)
end

get_mass_matrix(::Identity, ::Union{InPlace, OutOfPlace}, models, args...) = I

function get_mass_matrix(::Union{Constant, Varying}, ::InPlace, models)
    Nu = number_of_states(model)
    M = zeros(Nu,Nu)
    return get_mass_matrix!(models, M)
end

function get_mass_matrix(::Union{Constant, Varying}, ::InPlace, models, u, y, p, t)
    Nu = number_of_states(model)
    M = similar(u, Nu, Nu)
    return get_mass_matrix!(models, M, u, y, p, t)
end

function get_mass_matrix(::Constant, ::OutOfPlace, model::AbstractModel, u, y, p, t)
    return get_mass_matrix(model)
end

function get_mass_matrix(::Constant, ::OutOfPlace, models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    M = initialize_mass_matrix(models)

    # calculate input jacobian
    D = input_jacobian(models)

    # calculate input mass matrix
    My = get_input_mass_matrix(models...)

    return M + D*My
end

function get_mass_matrix(::Varying, ::OutOfPlace, models::NTuple{N, AbstractModel},
    u, y, p, t) where N

    # initialize mass matrix
    M = initialize_mass_matrix(models, u, y, p, t)

    # calculate input jacobian
    D = input_jacobian(models, u, y, p, t)

    # calculate input mass matrix
    My = get_input_mass_matrix(models..., u, p, t)

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

function input_jacobian(models::NTuple{N,AbstractModel}) where N
    static_input_jacobian(models)
end

function input_jacobian(models::NTuple{N,AbstractModel}, u, y, p, t) where N
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
        $(Di[1]) = vcat(get_input_jacobian(models[1]), $(Dij[1, 2:end]...))
    end
    for i = 2:N
        expr = quote
            $expr
            $(Di[i]) = vcat($(Dij[i, 1:i-1]...), get_input_jacobian(models[$i]),
                $(Dij[i, i+1:end]...))
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
    get_mass_matrix!(models, M)
    get_mass_matrix!(models, M, u, y, p, t)

In-place version of `get_mass_matrix`.
"""
get_mass_matrix!

function get_mass_matrix!(model::AbstractModel, M, u, y, p, t)
    return get_mass_matrix!(mass_matrix_state(model), inplaceness(model), model,
        M, u, y, p, t)
end

function get_mass_matrix!(models::NTuple{N,AbstractModel}, M, args...; kwargs...) where N
    return get_mass_matrix!(mass_matrix_state(models), inplaceness(models),
        models, M, args...; kwargs...)
end

function get_mass_matrix!(::Identity, ::Union{InPlace, OutOfPlace}, M,
    models::NTuple{N,AbstractModel}, args...; kwargs...) where N
    M .= 0
    for i = 1:N
        M[i,i] = 1
    end
    return M
end

function get_mass_matrix!(::Union{Constant, Varying}, ::OutOfPlace, models, M,
    args...; kwargs...)
    return M .= get_mass_matrix(models, args...)
end

function get_mass_matrix!(::Constant, ::InPlace, model::AbstractModel, M, u, y,
    p, t; kwargs...)
    return get_mass_matrix!(model, M; kwargs...)
end

function get_mass_matrix!(::Constant, ::InPlace, models::NTuple{N, AbstractModel},
    M; My = similar(M, number_of_inputs(models), number_of_states(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)

    # calculate input mass matrix
    get_input_mass_matrix!(models, My)

    # calculate mass matrix
    for i = 1:length(models)
        D = get_input_jacobian(models[i])
        for j = 1:length(models)
            Mij = view(M, iu[i], iu[j])
            Myij = view(My, iy[i], iu[j])
            if i == j
                get_mass_matrix!(models, Mij)
                if linear_input_dependence(models[i])
                    mul!(Mij, D, Myij, 1, 1)
                end
            else
                if linear_input_dependence(models[i])
                    mul!(Mij, D, Myij)
                else
                    Mij .= 0
                end
            end
        end
    end

    return M
end

function get_mass_matrix!(::Varying, ::InPlace, models::NTuple{N, AbstractModel},
    M; My = similar(M, number_of_inputs(models), number_of_states(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input mass matrix
    get_input_mass_matrix!(models, My, u, p, t)

    # calculate mass matrix
    for i = 1:length(models)
        ui = view(u, iu[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        D = get_input_jacobian(models[i], ui, yi, pi, t)
        for j = 1:length(models)
            Mij = view(M, iu[i], iu[j])
            Myij = view(My, iy[i], iu[j])
            if i == j
                get_mass_matrix!(models, Mij, ui, yi, pi, t)
                if linear_input_dependence(models[i])
                    mul!(Mij, D, Myij, 1, 1)
                end
            else
                if linear_input_dependence(models[i])
                    mul!(Mij, D, Myij)
                else
                    Mij .= 0
                end
            end
        end
    end

    return M
end

"""
    get_rates(models, u, y, p, t)

Calculate the (mass matrix multiplied) state rates for the specified model or
models.
"""
get_rates

function get_rates(model::AbstractModel, u, y, p, t)
    return get_rates(inplaceness(model), model, u, y, p, t)
end

function get_rates(models::NTuple{N,AbstractModel}, args...) where N
    return get_rates(inplaceness(models), models, args...)
end

function get_rates(::InPlace, models, u, y, p, t)
    du = similar(u)
    get_rates!(models, du, u, y, p, t)
    return du
end

function get_rates(::OutOfPlace, models::NTuple{N,AbstractModel}, u, y, p, t) where N
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
    get_rates!(models, du, u, y, p, t)

In-place version of [`get_rates`](@ref)
"""
get_rates!

function get_rates!(model::AbstractModel, du, u, y, p, t)
    return get_rates!(inplaceness(model), model, du, u, y, p, t)
end

function get_rates!(models::NTuple{N,AbstractModel}, du, u, y, p, t) where N
    return get_rates!(inplaceness(models), models, du, u, y, p, t)
end

function get_rates!(::OutOfPlace, models, du, u, y, p, t)
    return du .= get_rates(models, u, y, p, t)
end

function get_rates!(::InPlace, models::NTuple{N,AbstractModel}, du, u, y, p, t) where N

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    for i = 1:length(models)
        vdu = view(du, iu[i])
        vu = view(u, iu[i])
        vy = view(y, iy[i])
        vp = view(p, ip[i])
        get_rates!(models[i], vdu, vu, vy, vp, t)
    end

    return du
end

"""
    get_state_jacobian(models, u, y, p, t)

Calculate the jacobian with respect to the states for the specified model or models.
"""
get_state_jacobian

function get_state_jacobian(model::AbstractModel, u, y, p, t)
    if isinplace(model)
        N = number_of_states(model)
        J = zeros(N,N)
        get_state_jacobian!(model, J)
    else
        if has_state_jacobian(model)
            error("Required function `get_state_jacobian(model)` not defined "*
                "for model type $(typeof(model))")
        else
            N = number_of_states(model)
            J = SMatrix{N,N}(I)
        end
    end
    return J
end

function get_state_jacobian(models::NTuple{N,AbstractModel}, u, y, p, t) where N

    if isinplace(models)
        J = similar(u, (length(u), length(u)))
        get_state_jacobian!(models, J, u, y, p, t)
    else
        # get dimensions (these should be inferrable)
        Nu = number_of_states.(models)
        Ny = number_of_inputs.(models)
        Np = number_of_parameters.(models)

        # state variable indices
        iu2 = cumsum(Nu)
        iu1 = iu2 .- Nu .+ 1
        iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), length(models))

        # input variable indices
        iy2 = cumsum(Ny)
        iy1 = iy2 .- Ny .+ 1
        iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), length(models))

        # parameter indices
        ip2 = cumsum(Np)
        ip1 = ip2 .- Np .+ 1
        ip = ntuple(i->SVector{Np[i]}(ip1[i]:ip2[i]), length(models))

        # calculate input jacobian
        Jy = get_input_state_jacobian(models, u, p, t)

        # calculate jacobian
        Nutot = number_of_states(models)
        J = SMatrix{0,Nutot,Float64}()
        for i = 1:length(models)
            ui = u[iu[i]]
            yi = y[iy[i]]
            pi = p[ip[i]]
            D = get_input_jacobian(models[i], ui, yi, pi, t)
            Nui = number_of_states(models[i])
            Ji = SMatrix{Nui,0,Float64}()
            for j = 1:length(models)
                Jyij = Jy[iy[i], iu[j]]
                if i == j
                    Jij = get_state_jacobian(models[i], ui, yi, pi, t)
                    Jij = Jij + D*Jyij
                else
                    Jij = D*Jyij
                end
                Ji = hcat(Ji, Jij)
            end
            J = vcat(J, Ji)
        end
    end

    return J
end

"""
    get_state_jacobian!(models, J, u, y, p, t)

In-place version of [`get_state_jacobian`](@ref)
"""
get_state_jacobian!

function get_state_jacobian!(model::AbstractModel, J, u, y, p, t)
    if isinplace(model)
        error("Required method definition for `get_state_jacobian!` not defined for "*
            "model type $(typeof(model))")
    else
        J .= get_state_jacobian!(model, u, y, p, t)
    end
end

function get_state_jacobian!(models::NTuple{N,AbstractModel}, J, u, y, p, t;
    Jy = zeros(number_of_inputs(models), number_of_states(models))) where N

    if isinplace(models)
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
        get_input_state_jacobian!(models, Jy, u, p, t)

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
                    get_state_jacobian!(models, Jij, ui, yi, pi, t)
                    mul!(Jij, D, Jyij, 1, 1)
                else
                    mul!(Jij, D, Jyij)
                end
            end
        end
    else
        J .= get_state_jacobian(models, u, y, p, t)
    end

    return J
end

"""
    get_input_jacobian(models, J, u, y, p, t)

Calculate the jacobian with respect to the inputs for the specified model or models.
"""
function get_input_jacobian(model::AbstractModel, J, u, y, p, t)
    if isinplace(model)
        get_input_jacobian!(model, J, u, y, p, t)
    else
        error("Required method definition for `get_state_jacobian!` not defined for "*
            "model type $(typeof(model))")
    end
end

"""
    inplace_inputs(models)

Return a flag indicating whether the inputs corresponding to the model or
combination of models are calculated in-place.

This property should be inferrable from the model types.
"""
function inplace_inputs(models::NTuple{N,AbstractModel}) where N
    inplace_inputs(typeof.(models))
end

function inplace_inputs(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    inplace_inputs(T.parameters...)
end

"""
    has_input_mass_matrix(models)

Return a flag indicating whether the input function for the provided combination
of models has a mass matrix.

This property should be inferrable from the model types.
"""
function has_input_mass_matrix(models::NTuple{N,AbstractModel}) where N
    has_input_mass_matrix(typeof.(models))
end

function has_input_mass_matrix(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    has_input_mass_matrix(T.parameters...)
end

"""
    constant_input_mass_matrix(models)

Return a flag indicating whether the input function mass matrix is constant for
the provided combination of models

This property should be inferrable from the model type.
"""
function constant_input_mass_matrix(models::NTuple{N,AbstractModel}) where N
    constant_input_mass_matrix(typeof.(models))
end

function constant_input_mass_matrix(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    constant_input_mass_matrix(T.parameters...)
end

"""
    defined_input_state_jacobian(models)

Return a flag indicating whether the jacobian of the input function for the
provided combination of models is manually defined.
"""
function defined_input_state_jacobian(models::NTuple{N,AbstractModel}) where N
    defined_input_state_jacobian(typeof.(models))
end

function defined_input_state_jacobian(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    defined_input_state_jacobian(T.parameters...)
end

"""
    get_input_mass_matrix(models)
    get_input_mass_matrix(models, u, p, t)

Calculate the input function mass matrix for the specified combination of models.
"""
get_input_mass_matrix

function get_input_mass_matrix(models::NTuple{N,AbstractModel}) where N
    if inplace_inputs(models)
        Nu = number_of_states(models)
        Ny = number_of_inputs(models)
        My = zeros(Ny, Nu)
        get_input_mass_matrix!(models..., My)
    else
        My = get_input_mass_matrix(models...)
    end
    return My
end

function get_input_mass_matrix(models::NTuple{N,AbstractModel}, u, p, t) where N
    if inplace_inputs(models)
        Nu = number_of_states(models)
        Ny = number_of_inputs(models)
        My = similar(u, Ny, Nu)
        get_input_mass_matrix!(models..., My, u, p, t)
    else
        My = get_input_mass_matrix(models..., u, p, t)
    end
    return My
end

"""
    get_input_mass_matrix!(models, M)
    get_input_mass_matrix!(models, M, u, p, t)

In-place version of [`get_input_mass_matrix`](@ref).
"""
get_input_mass_matrix!

function get_input_mass_matrix!(models::NTuple{N,AbstractModel}, My) where N
    if inplace_inputs(models)
        get_input_mass_matrix!(models..., My)
    else
        My .= get_input_mass_matrix(models...)
    end
    return My
end

function get_input_mass_matrix!(models::NTuple{N,AbstractModel}, My, u, p, t) where N
    if inplace_inputs(models)
        get_input_mass_matrix!(models..., My, u, p, t)
    else
        My .= get_input_mass_matrix(models..., u, p, t)
    end
    return My
end

"""
    get_inputs(models, u, p, t)

Calculate the inputs to the specified combination of models.
"""
function get_inputs(models::NTuple{N,AbstractModel}, u, p, t) where N
    if inplace_inputs(models)
        Ny = number_of_inputs(models)
        y = similar(u, Ny)
        get_inputs!(models..., y, u, p, t)
    else
        y = get_inputs(models..., u, p, t)
    end
    return y
end

"""
    get_inputs!(models::NTuple{N,AbstractModel}, y, u, p, t) where N

In-place version of [`get_inputs`](@ref)
"""
function get_inputs!(models, y, u, p, t)
    if inplace_inputs
        get_inputs!(models..., y, u, p, t)
    else
        y .= get_inputs(models..., u, p, t)
    end
end

"""
    get_input_state_jacobian(models::NTuple{N,AbstractModel}, u, p, t) where N

Calculate the jacobian of the input function with respect to the state variables
for the specified combination of models.
"""
function get_input_state_jacobian(models::NTuple{N,AbstractModel}, u, p, t) where N
    if inplace_inputs(models)
        Nu = number_of_states(models)
        Ny = number_of_inputs(models)
        Jy = similar(u, Ny, Nu)
        get_input_state_jacobian!(models..., Jy, u, p, t)
    else
        Jy = get_input_mass_matrix(models..., u, p, t)
    end
    return Jy
end

"""
    get_input_state_jacobian!(models, J, u, p, t)

In-place version of [`get_input_state_jacobian`](@ref)
"""
function get_input_state_jacobian!(models::NTuple{N,AbstractModel}, Jy, u, p, t) where N
    if inplace_inputs(models)
        get_input_state_jacobian!(models..., Jy, u, p, t)
    else
        Jy .= get_input_state_jacobian(models..., u, p, t)
    end
    return Jy
end

# """
#     mass_matrix_operator!(models)
#
# Return the mass matrix corresponding to the specified combination of models as
# a DiffEqArrayOperator.
# """
# function mass_matrix_operator(models)
#
#     # problem dimensions
#     Nu = number_of_states(models)
#     Ny = number_of_inputs(models)
#
#     # initialize mass matrices
#     M = zeros(Nu, Nu)
#     My = zeros(Ny, Nu)
#
#     function update_func(M, u, p, t)
#         # update type of `My`
#         My = convert(typeof(M), My)
#         # calculate inputs
#
#          -> get_mass_matrix!(models, M, u, y, p, t; My = convert(typeof(M), My))
#     end
#
#     return DiffEqArrayOperator(M; update_func)
# end

# function mass_matrix_operator(models)
#
#     # initialize mass matrix
#     M = zeros(number_of_states(models),
#
#     # construct update function
#     if linear_input_dependence(stru)
#         Mas, Maa = init_load_mass_matrices(models)
#         update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t,
#             convert(typeof(M), Mas), convert(typeof(M), Maa))
#     else
#         update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t)
#     end
#
#     return DiffEqArrayOperator(M; update_func)
# end

# """
#     rate_function(models; input_cache = zeros(number_of_inputs(models)))
#
# Return a function that calculates the (mass matrix multiplied) state rates for
# the specified combination of models.
# """
# function rate_function(models; input_cache = zeros(number_of_inputs(models)))
#
#     # update inputs
#     if isinplace(models)
#
#
#
#     if isinplace(models)
#         f = (du, u, p, t) -> get_rates!(models, du, u, p, t)
#     else
#         f = (u, p, t) -> get_rates(models, u, p, t)
#     end
#
#     return f
# end

"""
    jacobian_function(models)

Return a function that calculates the jacobian for the specified combination of
models.
"""
function jacobian_function(models)

    if isinplace(models)
        jac = (J, u, p, t) -> get_state_jacobian!(models, J, u, p, t)
    else
        jac = (u, p, t) -> get_state_jacobian(models, u, p, t)
    end

    return jac
end

"""
    ODEFunction(models)

Construct an ODEFunction corresponding to the specified model or models which
may be solved using DifferentialEquations.
"""
function ODEFunction(models::NTuple{N,T}) where {N, T <: AbstractModel}

    # problem dimensions
    Nu = number_of_states(models)
    Ny = number_of_inputs(models)
    Np = number_of_parameters(models)

    # determine whether the problems
    iip = isinplace(models)

    # create cache variables
    ucache = fill(NaN, Nu)
    ycache = fill(NaN, Ny)
    pcache = fill(NaN, Np)
    tcache = fill(NaN)

    # construct in-place or out-of-place rate function
    if isinplace(models)
        f = (du, u, p, t) -> begin
            # check if we can use the cache variables (no custom types)
            if eltype(du) <: Float64
                # update the cache variables
                if (u != ucache) && (p != pcache) && (t != tcache[])
                    # calculate and store new model inputs
                    get_inputs!(models, ycache, u, p, t;
                        uperm = ucache, pperm = pcache)
                    # store current input arguments
                    ucache .= u
                    pcache .= p
                    tcache .= t
                end
                # use the cached model inputs
                y = ycache
            else
                # calculate model inputs (out-of-place to accomodate custom type)
                y = get_inputs(models, u, p, t)
            end
            # calculate mass matrix multiplied state rates
            get_rates!(models, du, u, y, p, t)
        end
    else
        f = (du, u, p, t) -> begin
            # check if we can use the cache variables (no custom types)
            if eltype(du) <: Float64
                # update the cache variables
                if (u != ucache) && (p != pcache) && (t != tcache[])
                    # calculate and store new model inputs
                    y = get_inputs(models, u, p, t)
                    # store current input arguments
                    ucache .= u
                    pcache .= p
                    tcache .= t
                end
                # use statically sized version of the cached model inputs
                y = SVector{Ny,Float64}(ycache)
            else
                # calculate model inputs (out-of-place to accomodate custom type)
                y = get_inputs(models, u, p, t)
            end
            # calculate mass matrix multiplied state rates
            du = get_rates(models, u, y, p, t)
        end
    end

    # construct mass matrix (or mass matrix operator)
    if has_mass_matrix(models)
        if constant_mass_matrix(models)
            mass_matrix = get_mass_matrix(models)
        else
            # initialize mass matrix
            M = zeros(Nu, Nu)
            # initialize cached input mass matrix
            My_cache = zeros(Ny, Nu)
            # construct update function
            update_func = (M, u, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if eltype(M) <: Float64
                    # update the cache variables
                    if (u != ucache) && (p != pcache) && (t != tcache[])
                        # calculate and store new model inputs
                        get_inputs!(models, ycache, u, p, t;
                            uperm = ucache, pperm = pcache)
                        # store current input arguments
                        ucache .= u
                        pcache .= p
                        tcache .= t
                    end
                    # use the cached model inputs
                    y = ycache
                else
                    # calculate model inputs (out-of-place to accomodate custom type)
                    y = get_inputs(models, u, p, t)
                end
                # update type of `My`
                My = convert(typeof(M), My_cache)
                # calculate inputs
                get_mass_matrix!(models, M, u, y, p, t; My = convert(typeof(M), My))
            end
            # construct mass matrix operator
            mass_matrix = DiffEqArrayOperator(M; update_func)
        end
    else
        mass_matrix = I
    end

    # construct jacobian function
    if defined_jacobian(models)
        if isinplace(models)
            jac = (J, u, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if eltype(du) <: Float64
                    # update the cache variables
                    if (u != ucache) && (p != pcache) && (t != tcache[])
                        # calculate and store new model inputs
                        get_inputs!(models, ycache, u, p, t;
                            uperm = ucache, pperm = pcache)
                        # store current input arguments
                        ucache .= u
                        pcache .= p
                        tcache .= t
                    end
                    # use the cached model inputs
                    y = ycache
                else
                    # calculate model inputs (out-of-place to accomodate custom type)
                    y = get_inputs(models, u, p, t)
                end
                # calculate jacobian
                get_state_jacobian!(models, J, u, y, p, t)
            end
        else
            jac = (u, p, t) -> begin
                # check if we can use the cache variables (no custom types)
                if eltype(du) <: Float64
                    # update the cache variables
                    if (u != ucache) && (p != pcache) && (t != tcache[])
                        # calculate and store new model inputs
                        y = get_inputs(models, u, p, t)
                        # store current input arguments
                        ucache .= u
                        pcache .= p
                        tcache .= t
                    end
                    # use statically sized version of the cached model inputs
                    y = SVector{Ny,Float64}(ycache)
                else
                    # calculate model inputs (out-of-place to accomodate custom type)
                    y = get_inputs(models, u, p, t)
                end
                # calculate jacobian
                J = get_state_jacobian(models, u, y, p, t)
            end
        end
    else
        jac = nothing # let DifferentialEquations construct the jacobian
    end

    # construct and return an ODEFunction
    return ODEFunction{iip}(f; mass_matrix, jac)
end
