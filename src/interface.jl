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

This property should be inferrable from the model types for in-place models.
"""
number_of_states

function number_of_states(model::AbstractModel)
    error("Required method definition for `number_of_states` not defined for "*
        "model type $(typeof(model))")
end

function number_of_states(models)
    N = 0
    for model in models
        N += number_of_states(model)
    end
    return N
end

"""
    number_of_inputs(models)

Return the total number of inputs corresponding to the model or models.

This property should be inferrable from the model types for in-place models.
"""
number_of_inputs

function number_of_inputs(model::AbstractModel)
    error("Required method definition for `number_of_inputs` not defined for "*
        "model type $(typeof(model))")
end

function number_of_inputs(models)
    N = 0
    for model in models
        N += number_of_inputs(model)
    end
    return N
end

"""
    isinplace(models)

Return a flag indicating whether the state rates corresponding to the model or
combination of models are calculated in-place.

This property should be inferrable from the model type.
"""
isinplace

function isinplace(model::AbstractModel)
    error("Required method definition for `isinplace` not defined for "*
        "model type $(typeof(model))")
end

function isinplace(models)
    inplace = inplace_input(models)
    for model in models
        inplace = inplace || isinplace(model)
    end
    return inplace
end

"""
    has_mass_matrix(models)

Return a flag indicating whether the governing differential equations for the
model or combination or models has a mass matrix.

This property should be inferrable from the model type.
"""
has_mass_matrix

function has_mass_matrix(model::AbstractModel)
    error("Required method definition for `has_mass_matrix` not defined for "*
        "model type $(typeof(model))")
end

function has_mass_matrix(models)
    result = has_input_mass_matrix(models)
    for model in models
        result = result || has_mass_matrix(model)
    end
    return result
end

"""
    constant_mass_matrix(models)

Return a flag indicating whether the mass matrix for the model or combination of
models is constant.

This property should be inferrable from the model type.
"""
constant_mass_matrix

function constant_mass_matrix(model::AbstractModel)
    if has_mass_matrix(model)
        error("Required method definition for `constant_mass_matrix` not defined "*
            "for model type $(typeof(model))")
    else
        return true
    end
end

function constant_mass_matrix(models)
    result = constant_input_mass_matrix(models)
    for model in models
        result = result && constant_input_jacobian(model)
        result = result && constant_mass_matrix(model)
    end
    return result
end

"""
    linear_input_dependence(models)

Return a flag indicating whether the state rates corresponding to the model or
combination of models has a linear dependence on the inputs
"""
linear_input_dependence

linear_input_dependence(model::AbstractModel) = false

function linear_input_dependence(models)
    result = true
    for model in models
        result = result && linear_input_dependence(model)
    end
    return result
end

"""
    defined_state_jacobian(models)

Return a flag indicating whether the jacobian of the state rates with respect
to the states is manually defined.
"""
defined_state_jacobian

defined_state_jacobian(model::AbstractModel) = false

function defined_state_jacobian(models)
    result = defined_input_state_jacobian(models)
    for model in models
        result = result && defined_state_jacobian(model)
        result = result && defined_input_jacobian(model)
    end
    return result
end

"""
    defined_input_jacobian(model)

Return a flag indicating whether the jacobian of the state rates with respect
to the inputs is defined. Defaults to false.
"""
defined_input_jacobian(model::AbstractModel) = false

"""
    constant_input_jacobian(model)

Return a flag indicating whether the jacobian of the state rates with respect
to the inputs is constant.  Defaults to false.
"""
constant_input_jacobian(model::AbstractModel) = false

"""
    get_mass_matrix(models)
    get_mass_matrix(models, u, y, p, t)

Calculate the mass matrix for a model or combination of models.
"""
get_mass_matrix

function get_mass_matrix(model::AbstractModel)
    if isinplace(model)
        N = number_of_states(model)
        M = zeros(N,N)
        get_mass_matrix!(model, M)
    else
        if has_mass_matrix(model)
            if constant_mass_matrix(model)
                error("Required function `get_mass_matrix(model)` not defined "*
                    "for model type $(typeof(model))")
            else
                error("State, input, parameter, and time arguments are required when "*
                    "calling `get_mass_matrix` for model type $(typeof(model))")
            end
        else
            N = number_of_states(model) # should be inferrable
            M = SMatrix{N,N}(I)
        end
    end
    return M
end

function get_mass_matrix(models)

    if isinplace(models)
        N = number_of_states(models)
        M = zeros(N, N)
        get_mass_matrix!(models, M)
    else
        if has_mass_matrix(models)
            if constant_mass_matrix(models)
                # get dimensions (these should be inferrable)
                Nu = number_of_states.(models)
                Ny = number_of_inputs.(models)

                # state variable indices
                iu2 = cumsum(Nu)
                iu1 = iu2 .- Nu
                iu = ntuple(i->SVector{Nu[i]}(range(iu1[i], iu2[i]), length(models))

                # input variable indices
                iy2 = cumsum(Ny)
                iy1 = iy2 .- Ny
                iy = ntuple(i->SVector{Ny[i]}(range(iy1[i], iy2[i]), length(models))

                # get input mass matrix
                My = get_input_mass_matrix(models)

                # calculate mass matrix
                N = number_of_states(models)
                M = SMatrix{0,N,Float64}()
                for i = 1:length(models)
                    D = get_input_jacobian(models[i])
                    Mi = SMatrix{Nu[i],0,Float64}()
                    for j = 1:length(models)
                        Myij = My[iy[i], iu[j]]
                        if i == j
                            Mij = get_mass_matrix(models[i])
                            if linear_input_dependence(models[i])
                                Mij = Mij + D*Myij
                            end
                        else
                            if linear_input_dependence(models[i])
                                Mij = D*Myij
                            else
                                Mij = SMatrix{Nu[i],Nu[j]}(ntuple(i->0, Nu*Nj))
                            end
                        end
                        Mi = hcat(Mi, Mij)
                    end
                    M = vcat(M, Mi)
                end
            else
                error("State, input, parameter, and time arguments are required when "*
                    "calling `get_mass_matrix` for model types $(typeof.(models)) ")
            end
        else
            N = number_of_states(models) # this should be inferrable
            M = SMatrix{N,N}(I)
        end
    end

    return M
end

function get_mass_matrix(models, u, y, p, t)

    if isinplace(models)
        N = number_of_states(models)
        M = zeros(N, N)
        get_mass_matrix!(models, M, u, y, p, t)
    else
        if has_mass_matrix(models)
            # get dimensions (these should be inferrable)
            Nu = number_of_states.(models)
            Ny = number_of_inputs.(models)
            Np = number_of_parameters.(models)

            # state variable indices
            iu2 = cumsum(Nu)
            iu1 = iu2 .- Nu
            iu = ntuple(i->SVector{Nu[i]}(range(iu1[i], iu2[i]), length(models))

            # input variable indices
            iy2 = cumsum(Ny)
            iy1 = iy2 .- Ny
            iy = ntuple(i->SVector{Ny[i]}(range(iy1[i], iy2[i]), length(models))

            # parameter indices
            ip2 = cumsum(Np)
            ip1 = ip2 .- Np
            ip = ntuple(i->SVector{Np[i]}(range(ip1[i], ip2[i]), length(models))

            # get input mass matrix
            My = get_input_mass_matrix(models, u, p, t)

            # calculate mass matrix
            N = number_of_states(models)
            M = SMatrix{0,N,Float64}()
            for i = 1:length(models)
                D = get_input_jacobian(models[i], u[iu[i]], p[ip[i]], t)
                Mi = SMatrix{Nu[i],0,Float64}()
                for j = 1:length(models)
                    Myij = My[iy[i], iu[j]]
                    if i == j
                        Mij = get_mass_matrix(models[i], u[iu[i]], y[iy[i]], p[ip[i]], t)
                        if linear_input_dependence(models[i])
                            Mij = Mij + D*Myij
                        end
                    else
                        if linear_input_dependence(models[i])
                            Mij = D*Myij
                        else
                            Mij = SMatrix{Nu[i],Nu[j]}(ntuple(i->0, Nu*Nj))
                        end
                    end
                    Mi = hcat(Mi, Mij)
                end
                M = vcat(M, Mi)
            end
        else
            N = number_of_states(models) # this should be inferrable
            M = SMatrix{N,N}(I)
        end
    end

    return M
end

"""
    get_mass_matrix!(models, M)
    get_mass_matrix!(models, M, u, y, p, t)

In-place version of `get_mass_matrix`.
"""
get_mass_matrix!

function get_mass_matrix!(model::AbstractModel, M)
    if isinplace(model)
        if has_mass_matrix(model)
            if constant_mass_matrix(model)
                error("Required function `get_mass_matrix!(model, M)` not defined "*
                    "for model type $(typeof(model))")
            else
                error("State, input, parameter, and time arguments are required when "*
                    "calling `get_mass_matrix!` for model type $(typeof(model))")
            end
        else
            N = number_of_states(model)
            M .= 0
            for i = 1:N
                M[i,i] = 1
            end
        end
    else
        M .= get_mass_matrix(model)
    end
    return M
end

function get_mass_matrix!(models, M;
    My = zeros(number_of_inputs(models), number_of_states(models)))

    if isinplace(models)
        if has_mass_matrix(models)
            if constant_mass_matrix(models)
                # get dimensions
                Nu = number_of_states.(models)
                Ny = number_of_inputs.(models)

                # state variable indices
                iu2 = cumsum(Nu)
                iu1 = iu2 .- Nu
                iu = range.(iu1, iu2)

                # input variable indices
                iy2 = cumsum(Ny)
                iy1 = iy2 .- Ny
                iy = range.(iy1, iy2)

                # get input mass matrix
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

            else
                error("State, input, parameter, and time arguments are required when "*
                    "calling `get_mass_matrix!` for model types $(typeof.(models)) ")
            end
        else
            N = number_of_states(models)
            M .= 0
            for i = 1:N
                M[i,i] = 1
            end
        end
    else
        M .= get_mass_matrix(models)
    end

    return M
end

function get_mass_matrix!(models, M, u, y, p, t;
    My = zeros(number_of_inputs(models), number_of_states(models)))

    if isinplace(models)
        if has_mass_matrix(models)
            # get dimensions
            Nu = number_of_states.(models)
            Ny = number_of_inputs.(models)
            Np = number_of_parameters.(models)

            # state variable indices
            iu2 = cumsum(Nu)
            iu1 = iu2 .- Nu
            iu = range.(iu1, iu2)

            # input variable indices
            iy2 = cumsum(Ny)
            iy1 = iy2 .- Ny
            iy = range.(iy1, iy2)

            # parameter indices
            ip2 = cumsum(Np)
            ip1 = ip2 .- Np
            ip = range.(ip1, ip2)

            # get input mass matrix
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
        else
            N = number_of_states(models)
            M .= 0
            for i = 1:N
                M[i,i] = 1
            end
        end
    else
        M .= get_mass_matrix(models, u, y, p, t)
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
    if isinplace(model)
        du = similar(u)
        get_rates!(model, du, u, y, p, t)
    else
        error("Required method definition for `get_rates` not defined for "*
            "model type $(typeof(model))")
    end
    return du
end

function get_rates(models, u, y, p, t)

    if isinplace(models)
        du = similar(u)
        get_rates!(models, du, u, y, p, t)
    else
        # get dimensions (these should be inferrable)
        Nu = number_of_states.(models)
        Ny = number_of_inputs.(models)
        Np = number_of_parameters.(models)

        # state variable indices
        iu2 = cumsum(Nu)
        iu1 = iu2 .- Nu
        iu = ntuple(i->SVector{Nu[i]}(range(iu1[i], iu2[i]), length(models))

        # input variable indices
        iy2 = cumsum(Ny)
        iy1 = iy2 .- Ny
        iy = ntuple(i->SVector{Ny[i]}(range(iy1[i], iy2[i]), length(models))

        # parameter indices
        ip2 = cumsum(Np)
        ip1 = ip2 .- Np
        ip = ntuple(i->SVector{Np[i]}(range(ip1[i], ip2[i]), length(models))

        # calculate state rates
        du = SVector{0,Float64}()
        for i = 1:length(models)
            du = vcat(du, get_rates(models[i], u[iu[i]], y[iy[i]], p[ip[i]], t))
        end
    end

    return du
end

"""
    get_rates!(models, du, u, y, p, t)

In-place version of [`get_rates`](@ref)
"""
get_rates!

function get_rates!(model::AbstractModel, du, u, y, p, t)
    if isinplace(model)
        error("Required method definition for `get_rates!` not defined for "*
            "model type $(typeof(model))")
    else
        du .= get_rates(model, u, y, p, t)
    end
    return du
end

function get_rates!(models, du, u, y, p, t)

    if isinplace(models)
        # get dimensions
        Nu = number_of_states.(models)
        Ny = number_of_inputs.(models)
        Np = number_of_parameters.(models)

        # state variable indices
        iu2 = cumsum(Nu)
        iu1 = iu2 .- Nu
        iu = range.(iu1, iu2)

        # input variable indices
        iy2 = cumsum(Ny)
        iy1 = iy2 .- Ny
        iy = range.(iy1, iy2)

        # parameter indices
        ip2 = cumsum(Np)
        ip1 = ip2 .- Np
        ip = range.(ip1, ip2)

        for i = 1:length(models)
            vdu = view(du, iu[i])
            vu = view(u, iu[i])
            vy = view(y, iy[i])
            vp = view(p, ip[i])
            get_rates!(models[i], vdu, vu, vy, vp, t)
        end
    else
        du .= get_rates(models, u, y, p, t)
    end

    return du
end

"""
    get_state_jacobian(models, u, y, p, t)

Calculate the jacobian with respect to the states for the specified model or models.
"""
get_state_jacobian

function get_state_jacobian(models, u, y, p, t)

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
        iu1 = iu2 .- Nu
        iu = ntuple(i->SVector{Nu[i]}(range(iu1[i], iu2[i]), length(models))

        # input variable indices
        iy2 = cumsum(Ny)
        iy1 = iy2 .- Ny
        iy = ntuple(i->SVector{Ny[i]}(range(iy1[i], iy2[i]), length(models))

        # parameter indices
        ip2 = cumsum(Np)
        ip1 = ip2 .- Np
        ip = ntuple(i->SVector{Np[i]}(range(ip1[i], ip2[i]), length(models))

        # calculate input jacobian
        Jy = get_input_jacobian(models, u, p, t)

        # calculate jacobian
        N = number_of_states(models)
        J = SMatrix{0,N,Float64}()
        for i = 1:length(models)
            ui = u[iu[i]]
            yi = y[iy[i]]
            pi = p[ip[i]]
            D = get_input_jacobian(models[i], ui, yi, pi, t)
            Ji = SMatrix{Nu[i],0,Float64}()
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
        J .= get_state_jacobian!(model, u, y, p, t)
    else
        error("Required method definition for `get_state_jacobian!` not defined for "*
            "model type $(typeof(model))")
    end
end

function get_state_jacobian!(models, J, u, y, p, t;
    Jy = zeros(number_of_inputs(models), number_of_states(models)))

    if isinplace(models)
        # get dimensions
        Nu = number_of_states.(models)
        Ny = number_of_inputs.(models)
        Np = number_of_parameters.(models)

        # state variable indices
        iu2 = cumsum(Nu)
        iu1 = iu2 .- Nu
        iu = range.(iu1, iu2)

        # input variable indices
        iy2 = cumsum(Ny)
        iy1 = iy2 .- Ny
        iy = range.(iy1, iy2)

        # parameter indices
        ip2 = cumsum(Np)
        ip1 = ip2 .- Np
        ip = range.(ip1, ip2)

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
get_input_jacobian

function get_input_jacobian(model::AbstractModel, J, u, y, p, t)
    if isinplace(model)
        get_input_jacobian!(model, J, u, y, p, t)
    else
        error("Required method definition for `get_state_jacobian!` not defined for "*
            "model type $(typeof(model))")
    end
end


"""
    mass_matrix_operator!(models)

Return the mass matrix corresponding to the specified combination of models as
a DiffEqArrayOperator.
"""
mass_matrix_operator(models)

function mass_matrix_operator(models::Tuple{AerodynamicModel, StructuralModel})
    return mass_matrix_operator(models[1], models[2])
end

function mass_matrix_operator(models::Tuple{StructuralModel, AerodynamicModel})
    return mass_matrix_operator(models[2], models[1])
end

function mass_matrix_operator(aero::AerodynamicModel, stru::StructuralModel)

    # initialize mass matrix
    M = init_mass_matrix(models)

    # construct update function
    if linear_load_dependence(stru)
        Mas, Maa = init_load_mass_matrices(models)
        update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t,
            convert(typeof(M), Mas), convert(typeof(M), Maa))
    else
        update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t)
    end

    return DiffEqArrayOperator(M; update_func)
end

"""
    jacobian_function!(models)

Return a function that calculates the jacobian for the specified combination of
models.
"""
jacobian_function(models)

function jacobian_function(models::Tuple{AerodynamicModel, StructuralModel})
    return jacobian_function(models[1], models[2])
end

function jacobian_function(models::Tuple{StructuralModel, AerodynamicModel})
    return jacobian_function(models[2], models[1])
end


function jacobian_function(aero, stru)
    if isinplace(models)
        jac = (J, u, p, t) -> get_state_jacobian!(models, J, u, p, t)
    else
        jac = (u, p, t) -> get_state_jacobian(models, u, p, t)
    end
end

function mass_matrix_operator(aero::AerodynamicModel, stru::StructuralModel)

    # initialize mass matrix
    M = init_mass_matrix(models)

    # construct update function
    if linear_load_dependence(stru)
        Mas, Maa = init_load_mass_matrices(models)
        update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t,
            convert(typeof(M), Mas), convert(typeof(M), Maa))
    else
        update_func = (M, u, p, t) -> update_mass_matrix!(models, M, u, p, t)
    end

    return DiffEqArrayOperator(M; update_func)
end

"""
    ODEFunction(models)

Construct an ODEFunction corresponding to the specified model or models which
may be solved using DifferentialEquations.
"""
function ODEFunction(models)

    # is the combined model in place?
    iip = isinplace(models)

    # construct in-place or out-of-place ODE
    if iip
        f = (du, u, p, t) -> get_rates!(models, du, u, p, t)
    else
        f = (u, p, t) -> get_rates(models, u, p, t)
    end

    # construct mass matrix
    if identity_mass_matrix(models)
        mass_matrix = I
    else
        if constant_mass_matrix(models)
            mass_matrix = init_mass_matrix(models)
        else
            mass_matrix = mass_matrix_operator(models)
        end
    end

    # construct jacobian function
    if defined_jacobian(models)
        jac = jacobian_function(models)
    else
        jac = nothing # let DifferentialEquations construct the jacobian
    end

    # construct and return an ODEFunction
    return ODEFunction{iip}(f; mass_matrix, jac)
end
