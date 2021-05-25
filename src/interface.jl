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

function number_of_states(models::NTuple{N,AbstractModel}) where N
    sum(number_of_states.(models))
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

function number_of_inputs(models::NTuple{N,AbstractModel}) where N
    sum(number_of_inputs.(models))
end

"""
    number_of_parameters(models)

Return the total number of parameters corresponding to the model or models.

This property should be inferrable from the model types for in-place models.
"""
number_of_parameters

function number_of_parameters(model::AbstractModel)
    error("Required method definition for `number_of_parameters` not defined for "*
        "model type $(typeof(model))")
end

function number_of_parameters(models::NTuple{N,AbstractModel}) where N
    sum(number_of_parameters.(models))
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

function isinplace(models::Tuple)
    inplace_inputs(models) || any(isinplace.(models))
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
    has_input_mass_matrix(models) || any(has_mass_matrix.(models))
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
    return constant_input_mass_matrix(models) &&
        all(constant_input_jacobian.(models)) &&
        all(constant_mass_matrix.(models))
end

"""
    linear_input_dependence(models)

Return a flag indicating whether the state rates corresponding to the model or
combination of models has a linear dependence on the inputs
"""
linear_input_dependence

linear_input_dependence(model::AbstractModel) = false

function linear_input_dependence(models)
    return all(linear_input_dependence.(models))
end

"""
    defined_state_jacobian(models)

Return a flag indicating whether the jacobian of the state rates with respect
to the states is manually defined.
"""
defined_state_jacobian

defined_state_jacobian(model::AbstractModel) = false

function defined_state_jacobian(models)
    return defined_input_state_jacobian(models) &&
        all(defined_state_jacobian.(models)) &&
        all(defined_input_jacobian.(models))
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
        Nu = number_of_states(model)
        M = zeros(Nu,Nu)
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
            Nu = number_of_states(model)
            M = SMatrix{Nu,Nu}(I)
        end
    end
    return M
end

function get_mass_matrix(model::AbstractModel, u, y, p, t)
    if isinplace(model)
        Nu = number_of_states(model)
        M = similar(u, Nu, Nu)
        get_mass_matrix!(model, M, u, y, p, t)
    else
        if has_mass_matrix(model)
            if constant_mass_matrix(model)
                M = get_mass_matrix(model)
            else
                error("Required function `get_mass_matrix(model)` not defined "*
                    "for model type $(typeof(model))")
            end
        else
            Nu = number_of_states(model)
            M = SMatrix{Nu,Nu}(I)
        end
    end
    return M
end

function get_mass_matrix(models::NTuple{N,AbstractModel}) where N

    if isinplace(models)
        Nu = number_of_states(models)
        M = zeros(Nu, Nu)
        get_mass_matrix!(models, M)
    else
        if has_mass_matrix(models)
            if constant_mass_matrix(models)
                # get dimensions (these should be inferrable)
                Nu = number_of_states.(models)
                Ny = number_of_inputs.(models)

                # state variable indices
                iu2 = cumsum(Nu)
                iu1 = iu2 .- Nu .+ 1
                iu = ntuple(i->SVector{Nu[i]}(iu1[i]:iu2[i]), length(models))

                # input variable indices
                iy2 = cumsum(Ny)
                iy1 = iy2 .- Ny .+ 1
                iy = ntuple(i->SVector{Ny[i]}(iy1[i]:iy2[i]), length(models))

                # get input mass matrix
                My = get_input_mass_matrix(models)

                # calculate mass matrix
                Nutot = number_of_states(models)
                M = SMatrix{0,Nutot,Float64}()
                for i = 1:length(models)
                    D = get_input_jacobian(models[i])
                    Nui = number_of_states(models[i])
                    Mi = SMatrix{Nui,0,Float64}()
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
                                Mij = zero(D*Myij)
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
            Nu = number_of_states(models) # this should be inferrable
            M = SMatrix{Nu,Nu}(I)
        end
    end

    return M
end

function get_mass_matrix(models::NTuple{N,AbstractModel}, u, y, p, t) where N

    if isinplace(models)
        Nu = number_of_states(models)
        M = zeros(Nu, Nu)
        get_mass_matrix!(models, M, u, y, p, t)
    else
        if has_mass_matrix(models)

            # get dimensions (these should be inferrable)
            Nu = number_of_states.(models)
            Ny = number_of_inputs.(models)
            Np = number_of_parameters.(models)

            # extract state variables for each model
            iu2 = cumsum(Nu)
            iu1 = iu2 .- Nu .+ 1
            iu = UnitRange.(iu1, iu2)
            us = getindex.(Ref(u), iu)

            # extract inputs for each model
            iy2 = cumsum(Ny)
            iy1 = iy2 .- Ny .+ 1
            iy = UnitRange.(iy1, iy2)
            ys = getindex.(Ref(y), iy)

            # extract parameters for each model
            ip2 = cumsum(Np)
            ip1 = ip2 .- Np .+ 1
            ip = UnitRange.(ip1, ip2)
            ps = getindex.(Ref(p), ip)

            # calculate mass matrix for each model
            Mii = get_mass_matrix.(models, us, ys, ps, t)

            # check for linear load dependence
            linear = linear_load_dependence.(models)

            # calculate mass matrix corresponding to the input function
            My = get_input_mass_matrix(models, u, p, t)

            #TODO: Finish this.

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

function get_mass_matrix!(model::AbstractModel, M, u, y, p, t)
    if isinplace(model)
        if has_mass_matrix(model)
            if constant_mass_matrix(model)
                get_mass_matrix!(model, M)
            else
                error("Required function `get_mass_matrix!(model, M)` not defined "*
                    "for model type $(typeof(model))")
            end
        else
            Nu = number_of_states(model)
            M .= 0
            for i = 1:Nu
                M[i,i] = 1
            end
        end
    else
        M .= get_mass_matrix(model, u, y, p, t)
    end
    return M
end

function get_mass_matrix!(models::NTuple{N,AbstractModel}, M;
    My = zeros(number_of_inputs(models), number_of_states(models))) where N

    if isinplace(models)
        if has_mass_matrix(models)
            if constant_mass_matrix(models)
                # get dimensions
                Nu = number_of_states.(models)
                Ny = number_of_inputs.(models)

                # state variable indices
                iu2 = cumsum(Nu)
                iu1 = iu2 .- Nu
                iu = UnitRange.(iu1, iu2)

                # input variable indices
                iy2 = cumsum(Ny)
                iy1 = iy2 .- Ny
                iy = UnitRange.(iy1, iy2)

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
            Nu = number_of_states(models)
            M .= 0
            for i = 1:Nu
                M[i,i] = 1
            end
        end
    else
        M .= get_mass_matrix(models)
    end

    return M
end

function get_mass_matrix!(models::NTuple{N,AbstractModel}, M, u, y, p, t;
    My = zeros(number_of_inputs(models), number_of_states(models))) where N

    if isinplace(models)
        if has_mass_matrix(models)
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
            Nu = number_of_states(models)
            M .= 0
            for i = 1:Nu
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

function get_rates(models::NTuple{N,AbstractModel}, u, y, p, t) where N

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

function get_rates!(models::NTuple{N,AbstractModel}, du, u, y, p, t) where N

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
    inplace_inputs(models...)
end

"""
    has_input_mass_matrix(models)

Return a flag indicating whether the input function for the provided combination
of models has a mass matrix.

This property should be inferrable from the model types.
"""
function has_input_mass_matrix(models::NTuple{N,AbstractModel}) where N
    has_input_mass_matrix(models...)
end

"""
    constant_input_mass_matrix(models)

Return a flag indicating whether the input function mass matrix is constant for
the provided combination of models

This property should be inferrable from the model type.
"""
function constant_input_mass_matrix(models::NTuple{N,AbstractModel}) where N
    constant_input_mass_matrix(models...)
end

"""
    defined_input_state_jacobian(models)

Return a flag indicating whether the jacobian of the input function for the
provided combination of models is manually defined.
"""
function defined_input_state_jacobian(models::NTuple{N,AbstractModel}) where N
    defined_input_state_jacobian(models...)
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
#     if linear_load_dependence(stru)
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
