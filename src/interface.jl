"""
    assemble_model(;
        aerodynamic_model = nothing, 
        structural_model = nothing, 
        dynamics_model = nothing,
        control_surface_model = nothing,
        controller_model = nothing,
        coupling_model = nothing)

Assemble a coupled model from a collection of models.
"""
function assemble_model(;
    aerodynamic_model = nothing, 
    structural_model = nothing, 
    dynamics_model = nothing,
    control_surface_model = nothing,
    controller_model = nothing,
    coupling_model = nothing)

    # construct tuple of provided models
    models = ()
    for model in (aerodynamic_model, structural_model, dynamics_model, 
        control_surface_model, controller_model)

        if !isnothing(model)
            models = (models..., model)
        end
    
    end

    # define submodels using provided models
    submodels = Submodel.(models)

    # use default coupling if a coupling model is not provided
    coupling = isnothing(coupling_model) ? Coupling(models...) : coupling_model

    # construct coupled model
    return CoupledModel(submodels, coupling)
end

"""
    assemble_states(model;
        aerodynamic_states = nothing, 
        structural_states = nothing, 
        dynamics_states = nothing,
        control_surface_states = nothing,
        controller_states = nothing)

Assemble the state vector for a coupled model using the specified state variables.  State 
variables are initialized with zeros.
"""
function assemble_states(model::CoupledModel; kwargs...)

    x = zeros(number_of_states(model))

    return assemble_states!(x, model; kwargs...)
end

"""
    assemble_states!(x, model;
        aerodynamic_states = nothing, 
        structural_states = nothing, 
        dynamics_states = nothing,
        control_surface_states = nothing,
        controller_states = nothing)

In-place version of [`assemble_states`](@ref).  Only replaces the values corresponding to 
the specified state variables.
"""
function assemble_states!(x0, model::CoupledModel;
    aerodynamic_states = nothing, 
    structural_states = nothing, 
    dynamics_states = nothing,
    control_surface_states = nothing,
    controller_states = nothing)

    # construct tuple of provided state variables
    coupled_states = ()
    for states in (aerodynamic_states, structural_states, dynamics_states,
        control_surface_states, controller_states)

        if !isnothing(states)
            coupled_states = (coupled_states..., states)
        end
    
    end

    # populate state vector
    for i = 1:length(states)
        set_states!(x0, model, i; coupled_states[i]...)
    end

    return x0
end

"""
    assemble_parameters(model;
        aerodynamic_parameters = nothing, 
        structural_parameters = nothing, 
        dynamics_parameters = nothing,
        control_surface_parameters = nothing,
        controller_parameters = nothing,
        additional_parameters = nothing)

Assemble the parameter vector for a coupled model using the specified parameter variables.
"""
function assemble_parameters(model::CoupledModel; kwargs...)

    x = zeros(number_of_parameters(model))

    return assemble_parameters!(x, model; kwargs...)
end

"""
    assemble_parameters!(p, model;
        aerodynamic_parameters = nothing, 
        structural_parameters = nothing, 
        dynamics_parameters = nothing,
        control_surface_parameters = nothing,
        controller_parameters = nothing,
        additional_parameters = nothing)

In-place version of [`assemble_parameters`](@ref). Only replaces the values corresponding 
to the specified parameters.
"""
function assemble_parameters!(p, model::CoupledModel;
    aerodynamic_parameters = nothing, 
    structural_parameters = nothing, 
    dynamics_parameters = nothing,
    control_surface_parameters = nothing,
    controller_parameters = nothing,
    additional_parameters = nothing)

    # define coupled parameters
    coupled_parameters = ()
    for parameters in (aerodynamic_parameters, structural_parameters, dynamics_parameters,
        control_surface_parameters, controller_parameters)

        if !isnothing(parameters)
            coupled_parameters = (coupled_parameters..., parameters)
        end
    
    end

    # populate parameter vector
    for i = 1:length(parameters)
        set_parameters!(p, model, i; coupled_parameters[i]...)
    end
    set_additional_parameters!(p, model; additional_parameters...)

    return p
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

"""
    ODEFunction(model::AerostructuralDynamics.CoupledModel, p=nothing)

Construct an `ODEFunction` for a coupled model which may be used with 
`DifferentialEquations`.     
"""
function SciMLBase.ODEFunction(model::CoupledModel, p=nothing)

    # always construct in-place models, always compile
    iip = true
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
        xcache = fill(NaN, Nx)
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
                if (x != xcache) && (p != pcache) && (t != tcache[])
                    xcache .= x
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

    f = (dx, x, p, t) -> begin
        get_residual!(dx, model, FillArrays.Zeros(x), x, fy(x, p, t), p, t)
        return dx .*= -1
    end

    return f
end

function ode_mass_matrix(model::CoupledModel, fy, p=nothing)

    submodels = model.submodels
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

        # state rate vector
        dx = FillArrays.Zeros(Nx)

        # initialize mass matrix
        M = zeros(Nx, Nx)

        # construct update function

        # cached variables
        drdy_cache = fill(NaN, Nx, Ny)
        dyddx_cache = fill(NaN, Ny, Nx)

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

        mass_matrix = DiffEqArrayOperator(M; update_func)
    end

    return mass_matrix 
end

function ode_time_gradient(model::CoupledModel, fy, p=nothing)

    submodels = model.submodels
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :tgrad))) && iszero(coupling.tgrad)
        # time gradient is a zero vector
        tgrad = (dT, x, p, t) -> dT .= 0
    elseif all(isinvariant.(getproperty.(submodels, :tgrad))) && isinvariant(coupling.tgrad)
        # time gradient is independent of all inputs
        dT0 = get_time_gradient(model)
        tgrad = (dT, x, p, t) -> dT .= dT0
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :tgrad))) && 
        isconstant(coupling.tgrad)
        # time gradient only depends on the parameter vector
        dT0 = get_rate_jacobian(model, p)
        tgrad = (dT, x, p, t) -> dT .= dT0
    else
        # time gradient defined as a function of state variables, parameters, and time 

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
    end

    return tgrad
end

function ode_state_jacobian(model::CoupledModel, fy, p=nothing)

    submodels = model.submodels
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :statejac))) && iszero(coupling.statejac)
        # state jacobian is a zero vector
        jac = (J, x, p, t) -> J .= 0
    elseif all(isinvariant.(getproperty.(submodels, :statejac))) && isinvariant(coupling.statejac)
        # state jacobian is independent of all inputs
        J0 = get_state_jacobian(model)
        jac = (J, x, p, t) -> J .= J0
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :statejac))) && 
        isconstant(coupling.statejac)
        # state jacobian only depends on the parameter vector
        J0 = get_state_jacobian(model, p)
        jac = (J, x, p, t) -> J .= J0
    else
        # state jacobian defined as a function of state variables, parameters, and time 

        # problem dimensions
        Nx = number_of_states(model)
        Ny = number_of_inputs(model)

        # state rate vector
        dx = FillArrays.Zeros(Nx)

        # construct state jacobian function

        # cached variables
        drdy_cache = fill(NaN, Nx, Ny)
        dyddx_cache = fill(NaN, Ny, Nx)

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
    end

    return jac
end

function ode_parameter_jacobian(model::CoupledModel, fy, p=nothing)

    submodels = model.submodels
    coupling = model.coupling

    if all(iszero.(getproperty.(submodels, :paramjac))) && isidentity(coupling.paramjac)
        # parameter jacobian is a zero matrix
        jac = (pJ, x, p, t) -> pJ .= 0
    elseif all(isinvariant.(getproperty.(submodels, :paramjac))) && isinvariant(coupling.paramjac)
        # parameter jacobian is independent of all inputs
        pJ0 = get_parameter_jacobian(model)
        paramjac = (pJ, x, p, t) -> pJ .= pJ0
    elseif !isnothing(p) && all(isconstant.(getproperty.(submodels, :paramjac))) && 
        isconstant(coupling.paramjac)
        # time gradient only depends on the parameter vector
        pJ0 = get_parameter_jacobian(model, p)
        paramjac = (pJ, x, p, t) -> pJ .= pJ0
    else
        # construct time gradient function

        # problem dimensions
        Nx = number_of_states(model)
        Ny = number_of_inputs(model)

        # cached variables
        drdy_cache = fill(NaN, Nx, Ny)
        dyddx_cache = fill(NaN, Ny, Nx)

        # state rate vector
        dx = FillArrays.Zeros(Nx)

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

    end

    return pjac
end

"""
    DAEFunction(model::AerostructuralDynamics.CoupledModel, p=nothing)

Construct an `DAEFunction` for a coupled model which may be used with 
`DifferentialEquations`.     
"""
function SciMLBase.DAEFunction(model::CoupledModel, p=nothing)

    # always construct in-place models, always compile
    iip = true
    compile = true

    # construct coupling function
    fy = dae_coupling_function(model)

    # construct rate function
    f = dae_residual_function(model, fy)

    # return ODE function
    return DAEFunction{iip, compile}(f)
end

function dae_coupling_function(model::CoupledModel)

    # check whether cached variables should be used
    if isinplace(model)
        # model is inplace so use cached variables if possible

        # problem dimensions
        Nx = number_of_states(model)
        Ny = number_of_inputs(model)
        Np = number_of_parameters(model)

        # cached variables
        dxcache = fill(NaN, Nx)
        xcache = fill(NaN, Nx)
        ycache = fill(NaN, Ny)
        pcache = fill(NaN, Np)
        tcache = fill(NaN)

        # coupling function
        fy = (dx, x, p, t) -> begin
            # check if we can use the cache variables (no custom types)
            if promote_type(eltype(dx), eltype(x), eltype(p), typeof(t)) <: Float64
                # update the cache variables (if necessary)
                if (dx != dxcache) && (x != xcache) && (p != pcache) && (t != tcache[])
                    dxcache .= dx
                    xcache .= x
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
        fy = (dx, x, p, t) -> get_coupling_inputs(model, dx, x, p, t)
    end

    # return coupling function
    return fy
end

function dae_residual_function(model::CoupledModel, fy)

    return (resid, dx, x, p, t) -> get_residual!(resid, model, dx, x, fy(dx, x, p, t), p, t)
end