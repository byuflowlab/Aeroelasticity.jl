"""
    number_of_states(models)

Return the total number of states corresponding to the model or models.
"""
number_of_states

# default to returning zero states
number_of_states(::Type{TM}) where TM = 0
number_of_states(model::TM) where TM = number_of_states(TM)

# models with no states have... no states
number_of_states(::Type{T}) where T<:NoStateModel = 0

# combined models have concatenated state variables
function number_of_states(models::NTuple{N,AbstractModel}) where N
    sum(number_of_states.(models))
end

# combined models have concatenated state variables
function number_of_states(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_states.(T.parameters))
end

"""
    number_of_inputs(models)

Return the total number of inputs corresponding to the model or models.
"""
number_of_inputs

# default to returning zero inputs
number_of_inputs(::Type{TM}) where TM = 0
number_of_inputs(model::TM) where TM = number_of_inputs(TM)

# models with no states have no inputs
number_of_inputs(::Type{T}) where T<:NoStateModel = 0

# combined models have concatenated inputs
function number_of_inputs(models::NTuple{N,AbstractModel}) where N
    sum(number_of_inputs.(models))
end

# combined models have concatenated inputs
function number_of_inputs(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_inputs.(T.parameters))
end

"""
    number_of_parameters(models)

Return the total number of parameters corresponding to the model or models.
"""
number_of_parameters

# default to returning zero parameters
number_of_parameters(::Vararg{Type,N}) where N = 0
number_of_parameters(models...) = number_of_parameters(typeof.(models)...)

# combined models have concatenated parameters
function number_of_parameters(models::NTuple{N,AbstractModel}) where N
    sum(number_of_parameters.(models)) + number_of_parameters(typeof.(models)...)
end

# combined models have concatenated parameters
function number_of_parameters(::Type{T}) where T<:NTuple{N,AbstractModel} where N
    sum(number_of_parameters.(T.parameters)) + number_of_parameters(T.parameters...)
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
    My = get_input_mass_matrix(models)

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
    My = get_input_mass_matrix(models, u, p, t)

    return M + D*My
end

# calculate mass matrix from mass matrix product
function _get_mass_matrix(::Undefined, ::OutOfPlace, models::T, args...) where T
    N = number_of_states(T)
    dui = SVector(ntuple(j -> ifelse(1 == j, 1, 0), N))
    M = get_mass_matrix_product(models, dui, args...)
    for i = 2:N
        dui = SVector(ntuple(j -> ifelse(i == j, 1, 0), N))
        out = get_mass_matrix_product(models, dui, args...)
        M = vcat(M, out)
    end
    return M
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
    get_input_mass_matrix!(My, models)

    # calculate mass matrix
    for i = 1:length(models)
        D = get_input_jacobian(models[i])
        for j = 1:length(models)
            Mij = view(M, iu[i], iu[j])
            Myij = view(My, iy[i], iu[j])
            if i == j
                get_mass_matrix!(Mij, models)
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

# calculate mass matrix for a combination of models
function _get_mass_matrix!(M, ::Linear, ::InPlace,
    models::NTuple{N, AbstractModel}; My = similar(M, number_of_inputs(models),
    number_of_states(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input mass matrix
    get_input_mass_matrix!(My, models, u, p, t)

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
                get_mass_matrix!(Mij, models, ui, yi, pi, t)
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

# calculate mass matrix from mass matrix product
function _get_mass_matrix!(M, ::Undefined, ::OutOfPlace, models::T, args...) where T
    N = number_of_states(T)
    dui = SVector(ntuple(j -> ifelse(1 == j, 1, 0), N))
    get_mass_matrix_product!(view(M, i, :), models, dui, args...)
    for i = 2:N
        dui = SVector(ntuple(j -> ifelse(i == j, 1, 0), N))
        get_mass_matrix_product!(view(M, i, :), models, dui, args...)
    end
    return M
end

"""
    get_mass_matrix_product(models, du, u, y, p, t)

Calculate the mass matrix product with the state rates for the specified model
or models.
"""
function get_mass_matrix_product(models::T, du, u, y, p, t) where T
    return _get_mass_matrix_product(mass_matrix_type(T), inplaceness(T), models,
        du, u, y, p, t)
end

# dispatch to an in-place function
function _get_mass_matrix_product(::Any, ::InPlace, models, du)
    out = similar(du)
    get_mass_matrix_product!(out, models, du)
    return du
end

# dispatch to an in-place function
function _get_mass_matrix_product(::Any, ::InPlace, models, du, u, y, p, t)
    TF = promote_type(eltype(du), eltype(u), eltype(y), eltype(p), typeof(t))
    out = similar(du, TF)
    get_mass_matrix_product!(out, models, du, u, y, p, t)
    return out
end

# return an empty vector
function _get_mass_matrix_product(::Empty, ::OutOfPlace, models, args...)
    return SVector{0,Float64}()
end

# return the state rates
function _get_mass_matrix_product(::Identity, ::OutOfPlace, models, du, args...)
    return du
end

# dispatch to `get_mass_matrix_product` without arguments
function _get_mass_matrix_product(::Constant, ::OutOfPlace, models, du, u, y, p, t)
    return get_mass_matrix_product(models, du)
end

# calculate mass matrix for a combination of models
function _get_mass_matrix_product(::Constant, ::OutOfPlace, models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix product
    Mdu = initialize_mass_matrix_product(models, du)

    # calculate input jacobian
    D = input_jacobian(models)

    # calculate input mass matrix product
    Mydu = get_input_mass_matrix_product(models, du)

    return Mdu + D*Mydu
end

# calculate mass matrix for a combination of models
function _get_mass_matrix_product(::Linear, ::OutOfPlace, models::NTuple{N, AbstractModel},
    u, y, p, t) where N

    # initialize mass matrix product
    Mdu = initialize_mass_matrix_product(models, du, u, y, p, t)

    # calculate input jacobian
    D = get_input_jacobian(models, u, y, p, t)

    # calculate input mass matrix product
    Mydu = get_input_mass_matrix_product(models, du, u, p, t)

    return Mdu + D*Mydu
end

function initialize_mass_matrix_product(models::NTuple{N,AbstractModel}) where N
    initialize_static_mass_matrix_product(models)
end

function initialize_mass_matrix_product(models::NTuple{N,AbstractModel}, u, y, p, t) where N
    initialize_varying_mass_matrix_product(models, u, y, p, t)
end

@generated function initialize_static_mass_matrix_product(models::NTuple{N,AbstractModel}) where N

    # get indices of state, input, and parameter vectors for each model
    Nu = number_of_states.(models.parameters)
    iu2 = cumsum(Nu)
    iu1 = iu2 .- Nu .+ 1
    iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), N)

    # initialize state variable names
    outi = [Symbol("out", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(outi[i]) = get_mass_matrix_product(models[$i], du[$(iu[i])])
        end
    end

    expr = quote
        $expr
        out = vcat($(dui...))
    end

    return expr
end

@generated function initialize_varying_mass_matrix_product(models::NTuple{N,AbstractModel},
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

    # initialize state variable names
    outi = [Symbol("out", i) for i = 1:N]

    expr = :()
    for i = 1:N
        expr = quote
            $expr
            $(outi[i]) = get_mass_matrix_product(models[$i], du[$(iu[i])],
                u[$(iu[i])], y[$(iy[i])], p[$(ip[i])], t)
        end
    end

    expr = quote
        $expr
        out = vcat($(outi...))
    end

    return expr
end

"""
    get_mass_matrix_product!(out, models, du, u, y, p, t)

In-place version of [`get_rates`](@ref)
"""
function get_mass_matrix_product!(out, models::T, u, y, p, t) where T
    return _get_mass_matrix_product!(out, inplaceness(T), models, du, u, y, p, t)
end

# dispatch to an out-of-place function
function _get_mass_matrix_product!(out, ::Any, ::OutOfPlace, models, args...)
    return out .= get_mass_matrix_product(models, args...)
end

# return an empty vector
function _get_mass_matrix_product!(out, ::Empty, ::InPlace, models, args...)
    return out
end

# return the state rates
function _get_mass_matrix_product!(out, ::Identity, ::InPlace, models, du, args...)
    copyto!(out, du)
    return out
end

# dispatch to `get_mass_matrix_product` without arguments
function _get_mass_matrix_product!(out, ::Constant, ::InPlace, models, du, u, y, p, t)
    return get_mass_matrix_product!(out, models, du)
end

# calculate mass matrix for a combination of models
function _get_mass_matrix_product!(Mdu, ::Constant, ::InPlace,
    models::NTuple{N, AbstractModel}, du;
    Mydu = similar(Mdu, number_of_inputs(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)

    # calculate input mass matrix
    get_input_mass_matrix_product!(Mydu, models, du)

    # calculate mass matrix
    for i = 1:length(models)
        Mdui = view(Mdu, iu[i])
        Mydui = view(Mydu, iy[i])
        D = get_input_jacobian(models[i])

        # calculate mass matrix product contributions from self
        get_mass_matrix_product!(Mdui, models[i], view(u, iu[i]))

        # add mass matrix product contribution from other models
        mul!(Mdui, D, Mydui, 1, 1)
    end

    return Mdu
end

# calculate mass matrix product for a combination of models
function _get_mass_matrix_product!(Mdu, ::Linear, ::InPlace,
    models::NTuple{N, AbstractModel}; Mydu = similar(Mdu,
    number_of_inputs(models))) where N

    # get state and parameter indices
    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    # calculate input mass matrix
    get_input_mass_matrix_product!(Mydu, models, u, p, t)

    # calculate mass matrix
    for i = 1:length(models)
        dui = view(du, iu[i])
        ui = view(u, iu[i])
        yi = view(y, iy[i])
        pi = view(p, ip[i])
        Mdui = view(Mdu, iu[i])
        Mydui = view(Mydu, iy[i])
        D = get_input_jacobian(models[i], ui, yi, pi, t)

        # calculate mass matrix product contributions from self
        Mdui = view(out, iu[i])
        get_mass_matrix_product!(Mdui, models[i], dui, ui, yi, pi, t)

        # add mass matrix product contributions from other models
        D = get_input_jacobian(models[i], ui, yi, pi, t)
        Mydui = view(Mydu, iy[i])
        mul!(Mdui, D, Mydui, 1, 1)
    end

    return Mdu
end

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

# return the identity matrix
function _get_state_jacobian(::Identity, ::OutOfPlace, models::TM, args...) where TM
    Nu = number_of_states(TM)
    return SMatrix{Nu, Nu, Float64}(I)
end

# dispatch to the `get_state_jacobian` without arguments
function _get_state_jacobian(::Constant, ::OutOfPlace, models, u, y, p, t)
    return get_state_jacobian(models)
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}) where N

    # initialize mass matrix
    J = initialize_state_jacobian(models)

    # calculate input jacobian
    D = input_jacobian(models)

    # calculate input jacobian
    Jy = get_input_state_jacobian(models)

    return J + D*Jy
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace, models,
    u, y, p, t)

    return ForwardDiff.jacobian(u->get_rates(models, u, y, p, t), u)
end

# calculate state jacobian for a combination of models
function _get_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N, AbstractModel}, u, y, p, t) where N

    # initialize mass matrix
    J = initialize_state_jacobian(models, u, y, p, t)

    # calculate input jacobian
    D = get_input_jacobian(models, u, y, p, t)

    # calculate input jacobian
    Jy = get_input_state_jacobian(models, u, p, t)

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
    get_input_state_jacobian!(Jy, models)

    # calculate jacobian
    for i = 1:length(models)
        D = get_input_jacobian(models[i])
        for j = 1:length(models)
            Jij = view(J, iu[i], iu[j])
            Jyij = view(Jy, iy[i], iu[j])
            if i == j
                get_state_jacobian!(Jij, models)
                mul!(Jij, D, Jyij, 1, 1)
            else
                mul!(Jij, D, Jyij)
            end
        end
    end

    return J
end

# use automatic differentiation since a custom definition is absent
function _get_state_jacobian!(J, ::Union{Linear, Nonlinear}, ::InPlace,
    model::AbstractModel, u, y, p, t)

    return ForwardDiff.jacobian!(J, u->get_rates(models, u, y, p, t), u)
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
    get_input_state_jacobian!(Jy, models, u, p, t)

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
                get_state_jacobian!(Jij, models, ui, yi, pi, t)
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

# dispatch to the `get_input_jacobian` function without arguments
function _get_input_jacobian(::Constant, models, u, y, p, t)
    return get_input_jacobian(models)
end

# calculate input jacobian for a combination of models
function _get_input_jacobian(::Constant, models::NTuple{N,AbstractModel}) where N
    static_input_jacobian(models)
end

# use automatic differentiation since a custom definition is absent
function _get_input_jacobian(::Union{Linear, Nonlinear}, model::AbstractModel, u, y, p, t)

    return ForwardDiff.jacobian(y->get_rates(model, u, y, p, t), y)
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
    get_input_mass_matrix(models)
    get_input_mass_matrix(models, u, p, t)

Calculate the input function mass matrix for the specified combination of models.
"""
get_input_mass_matrix

function get_input_mass_matrix(models::TM; kwargs...) where TM
    return _get_input_mass_matrix(mass_matrix_type(TM.parameters...),
        inplaceness(TM.parameters...), models; kwargs...)
end

function get_input_mass_matrix(models::TM, u, p, t; kwargs...) where TM
    return _get_input_mass_matrix(mass_matrix_type(TM.parameters...),
        inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an in-place function
function _get_input_mass_matrix(::Any, ::InPlace, models; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = zeros(Ny,Nu)
    return get_input_mass_matrix!(M, models; kwargs...)
end

# dispatch to an in-place function
function _get_input_mass_matrix(::Any, ::InPlace, models, u, p, t; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = similar(u, Ny, Nu)
    return get_input_mass_matrix!(M, models, u, p, t; kwargs...)
end

# return an empty matrix
function _get_input_mass_matrix(::Empty, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)

    return SMatrix{Ny, 0, Float64}()
end

# return a matrix of zeros
function _get_input_mass_matrix(::Zeros, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)
    Nu = number_of_states(TM)

    return zeros(SMatrix{Ny, Nu, Float64})
end

# dispatch to the `get_input_mass_matrix` function without arguments
function _get_input_mass_matrix(::Constant, ::OutOfPlace, models, u, p, t; kwargs...)
    return get_input_mass_matrix(models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_mass_matrix(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    return get_input_mass_matrix(models...; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_mass_matrix(::Linear, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    return get_input_mass_matrix(models..., u, p, t; kwargs...)
end

"""
    get_input_mass_matrix!(My, models)
    get_input_mass_matrix!(My, models, u, p, t)

In-place version of [`get_input_mass_matrix`](@ref).
"""
get_input_mass_matrix!

function get_input_mass_matrix!(My, models::TM; kwargs...) where TM
    return _get_input_mass_matrix!(My, mass_matrix_type(TM.parameters...),
        inplaceness(TM.parameters...), models; kwargs...)
end

function get_input_mass_matrix!(My, models::TM, u, p, t; kwargs...) where TM
    return _get_input_mass_matrix!(My, mass_matrix_type(TM.parameters...),
        inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_input_mass_matrix!(My, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return My .= get_input_mass_matrix(models, args...; kwargs...)
end

# return an empty matrix
function _get_input_mass_matrix!(My, ::Empty, ::InPlace, models, args...; kwargs...)
    return My
end

# return a matrix of zeros
function _get_input_mass_matrix!(My, ::Zeros, ::InPlace, models, args...; kwargs...)
    return My .= 0
end

# dispatch to `get_input_mass_matrix!` without arguments
function _get_input_mass_matrix!(My, ::Constant, ::InPlace, models, u, p, t; kwargs...)
    return get_input_mass_matrix!(My, models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_mass_matrix!(My, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    get_input_mass_matrix!(My, models...; kwargs...)

    return My
end

# dispatch to the user-provided function for the specific combination of models
function get_input_mass_matrix!(My, ::Linear, ::InPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    get_input_mass_matrix!(My, models..., u, p, t; kwargs...)

    return My
end

"""
    get_inputs(models, u, p, t)

Calculate the inputs to the specified combination of models.
"""
function get_inputs(models::TM, u, p, t) where TM
    return _get_inputs(inplaceness(TM.parameters...), models, u, p, t)
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
    return _get_inputs!(y, inplaceness(T.parameters...), models, u, p, t)
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
    get_input_state_jacobian(models::NTuple{N,AbstractModel}, u, p, t) where N

Calculate the jacobian of the input function with respect to the state variables
for the specified combination of models.
"""
get_input_state_jacobian

function get_input_state_jacobian(models::TM; kwargs...) where TM
    return _get_input_state_jacobian(state_jacobian_type(TM.parameters...),
        inplaceness(TM.parameters...), models; kwargs...)
end

function get_input_state_jacobian(models::TM, u, p, t; kwargs...) where TM
    return _get_input_state_jacobian(state_jacobian_type(TM.parameters...),
        inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an in-place function
function _get_input_state_jacobian(::Any, ::InPlace, models; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = zeros(Ny,Nu)
    return get_input_state_jacobian!(M, models; kwargs...)
end

# dispatch to an in-place function
function _get_input_state_jacobian(::Any, ::InPlace, models, u, p, t; kwargs...)
    Ny = number_of_inputs(models)
    Nu = number_of_states(models)
    M = similar(u, Ny, Nu)
    return get_input_state_jacobian!(M, models, u, p, t; kwargs...)
end

# return an empty matrix
function _get_input_state_jacobian(::Empty, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)

    return SMatrix{Ny, 0, Float64}()
end

# return a matrix of zeros
function _get_input_state_jacobian(::Zeros, ::OutOfPlace, models::TM, args...;
    kwargs...) where TM

    Ny = number_of_inputs(TM)
    Nu = number_of_states(TM)

    return zeros(SMatrix{Ny, Nu, Float64})
end

# dispatch to the `get_input_state_jacobian` function without arguments
function _get_input_state_jacobian(::Constant, ::OutOfPlace, models, u, p, t; kwargs...)
    return get_input_state_jacobian(models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_state_jacobian(::Constant, ::OutOfPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    return get_input_state_jacobian(models...; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_state_jacobian(::Union{Linear, Nonlinear}, ::OutOfPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    return get_input_state_jacobian(models..., u, p, t; kwargs...)
end

"""
    get_input_state_jacobian!(J, models, u, p, t)

In-place version of [`get_input_state_jacobian`](@ref)
"""
get_input_state_jacobian!

function get_input_state_jacobian!(Jy, models::TM; kwargs...) where TM
    return _get_input_state_jacobian!(Jy, state_jacobian_type(TM.parameters...),
        inplaceness(TM.parameters...), models; kwargs...)
end

function get_input_state_jacobian!(Jy, models::TM, u, p, t; kwargs...) where TM
    return _get_input_state_jacobian!(Jy, state_jacobian_type(TM.parameters...),
        inplaceness(TM.parameters...), models, u, p, t; kwargs...)
end

# dispatch to an out-of-place function
function _get_input_state_jacobian!(Jy, ::Any, ::OutOfPlace, models, args...; kwargs...)
    return Jy .= get_input_state_jacobian(models, args...; kwargs...)
end

# return an empty matrix
function _get_input_state_jacobian!(Jy, ::Empty, ::InPlace, models, args...; kwargs...)
    return Jy
end

# return a matrix of zeros
function _get_input_state_jacobian!(Jy, ::Zeros, ::InPlace, models, args...; kwargs...)
    return Jy .= 0
end

# dispatch to `get_input_state_jacobian!` without arguments
function _get_input_state_jacobian!(Jy, ::Constant, ::InPlace, models, u, p, t; kwargs...)
    return get_input_state_jacobian!(Jy, models; kwargs...)
end

# dispatch to the user-provided function for the specific combination of models
function _get_input_state_jacobian!(Jy, ::Constant, ::InPlace,
    models::NTuple{N,AbstractModel}; kwargs...) where N

    get_input_state_jacobian!(Jy, models...; kwargs...)

    return Jy
end

# dispatch to the user-provided function for the specific combination of models
function get_input_state_jacobian!(Jy, ::Union{Linear, Nonlinear}, ::InPlace,
    models::NTuple{N,AbstractModel}, u, p, t; kwargs...) where N

    get_input_state_jacobian!(Jy, models..., u, p, t; kwargs...)

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
#
# """
#     jacobian_function(models)
#
# Return a function that calculates the jacobian for the specified combination of
# models.
# """
# function jacobian_function(models)
#
#     if isinplace(models)
#         jac = (J, u, p, t) -> get_state_jacobian!(J, models, u, p, t)
#     else
#         jac = (u, p, t) -> get_state_jacobian(models, u, p, t)
#     end
#
#     return jac
# end
#
# """
#     ODEFunction(models)
#
# Construct an ODEFunction corresponding to the specified model or models which
# may be solved using DifferentialEquations.
# """
# function ODEFunction(models::NTuple{N,T}) where {N, T <: AbstractModel}
#
#     # problem dimensions
#     Nu = number_of_states(models)
#     Ny = number_of_inputs(models)
#     Np = number_of_parameters(models)
#
#     # determine whether the problems
#     iip = isinplace(models)
#
#     # create cache variables
#     ucache = fill(NaN, Nu)
#     ycache = fill(NaN, Ny)
#     pcache = fill(NaN, Np)
#     tcache = fill(NaN)
#
#     # construct in-place or out-of-place rate function
#     if isinplace(models)
#         f = (du, u, p, t) -> begin
#             # check if we can use the cache variables (no custom types)
#             if eltype(du) <: Float64
#                 # update the cache variables
#                 if (u != ucache) && (p != pcache) && (t != tcache[])
#                     # calculate and store new model inputs
#                     get_inputs!(ycache, models, u, p, t;
#                         uperm = ucache, pperm = pcache)
#                     # store current input arguments
#                     ucache .= u
#                     pcache .= p
#                     tcache .= t
#                 end
#                 # use the cached model inputs
#                 y = ycache
#             else
#                 # calculate model inputs (out-of-place to accomodate custom type)
#                 y = get_inputs(models, u, p, t)
#             end
#             # calculate mass matrix multiplied state rates
#             get_rates!(du, models, u, y, p, t)
#         end
#     else
#         f = (du, u, p, t) -> begin
#             # check if we can use the cache variables (no custom types)
#             if eltype(du) <: Float64
#                 # update the cache variables
#                 if (u != ucache) && (p != pcache) && (t != tcache[])
#                     # calculate and store new model inputs
#                     y = get_inputs(models, u, p, t)
#                     # store current input arguments
#                     ucache .= u
#                     pcache .= p
#                     tcache .= t
#                 end
#                 # use statically sized version of the cached model inputs
#                 y = SVector{Ny,Float64}(ycache)
#             else
#                 # calculate model inputs (out-of-place to accomodate custom type)
#                 y = get_inputs(models, u, p, t)
#             end
#             # calculate mass matrix multiplied state rates
#             du = get_rates(models, u, y, p, t)
#         end
#     end
#
#     # construct mass matrix (or mass matrix operator)
#     if has_mass_matrix(models)
#         if constant_mass_matrix(models)
#             mass_matrix = get_mass_matrix(models)
#         else
#             # initialize mass matrix
#             M = zeros(Nu, Nu)
#             # initialize cached input mass matrix
#             My_cache = zeros(Ny, Nu)
#             # construct update function
#             update_func = (M, u, p, t) -> begin
#                 # check if we can use the cache variables (no custom types)
#                 if eltype(M) <: Float64
#                     # update the cache variables
#                     if (u != ucache) && (p != pcache) && (t != tcache[])
#                         # calculate and store new model inputs
#                         get_inputs!(ycache, models, u, p, t;
#                             uperm = ucache, pperm = pcache)
#                         # store current input arguments
#                         ucache .= u
#                         pcache .= p
#                         tcache .= t
#                     end
#                     # use the cached model inputs
#                     y = ycache
#                 else
#                     # calculate model inputs (out-of-place to accomodate custom type)
#                     y = get_inputs(models, u, p, t)
#                 end
#                 # update type of `My`
#                 My = convert(typeof(M), My_cache)
#                 # calculate inputs
#                 get_mass_matrix!(M, models, u, y, p, t; My = convert(typeof(M), My))
#             end
#             # construct mass matrix operator
#             mass_matrix = DiffEqArrayOperator(M; update_func)
#         end
#     else
#         mass_matrix = I
#     end
#
#     # construct jacobian function
#     if defined_jacobian(models)
#         if isinplace(models)
#             jac = (J, u, p, t) -> begin
#                 # check if we can use the cache variables (no custom types)
#                 if eltype(du) <: Float64
#                     # update the cache variables
#                     if (u != ucache) && (p != pcache) && (t != tcache[])
#                         # calculate and store new model inputs
#                         get_inputs!(ycache, models, u, p, t;
#                             uperm = ucache, pperm = pcache)
#                         # store current input arguments
#                         ucache .= u
#                         pcache .= p
#                         tcache .= t
#                     end
#                     # use the cached model inputs
#                     y = ycache
#                 else
#                     # calculate model inputs (out-of-place to accomodate custom type)
#                     y = get_inputs(models, u, p, t)
#                 end
#                 # calculate jacobian
#                 get_state_jacobian!(J, models, u, y, p, t)
#             end
#         else
#             jac = (u, p, t) -> begin
#                 # check if we can use the cache variables (no custom types)
#                 if eltype(du) <: Float64
#                     # update the cache variables
#                     if (u != ucache) && (p != pcache) && (t != tcache[])
#                         # calculate and store new model inputs
#                         y = get_inputs(models, u, p, t)
#                         # store current input arguments
#                         ucache .= u
#                         pcache .= p
#                         tcache .= t
#                     end
#                     # use statically sized version of the cached model inputs
#                     y = SVector{Ny,Float64}(ycache)
#                 else
#                     # calculate model inputs (out-of-place to accomodate custom type)
#                     y = get_inputs(models, u, p, t)
#                 end
#                 # calculate jacobian
#                 J = get_state_jacobian(models, u, y, p, t)
#             end
#         end
#     else
#         jac = nothing # let DifferentialEquations construct the jacobian
#     end
#
#     # construct and return an ODEFunction
#     return ODEFunction{iip}(f; mass_matrix, jac)
# end
