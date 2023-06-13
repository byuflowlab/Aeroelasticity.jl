"""
    CoupledModel{TR, TY, TP, TI, TS, TC, TM}

Base type for coupled models in Aeroelasticity.jl.
"""
struct CoupledModel{TR, TY, TP, TI, TS, TC, TM}
    fresid::TR
    finput::TY
    fparam::TP
    indices::TI
    rate_jacobian_sparsity::TS
    state_jacobian_sparsity::TS
    rate_jacobian_colorvec::TC
    state_jacobian_colorvec::TC
    identity_mass_matrix::Bool
    constant_mass_matrix::Bool
    mass_matrix::TM
end

"""
    CoupledModel(submodels, parameters, nstate=number_of_states.(submodels); kwargs...)

Define a coupled model from a collection of submodels.  Each submodel is defined as an
in-place implicit differential equation (i.e. `submodels[i](residᵢ, dxᵢ, xᵢ, yᵢ, pᵢ, t)`).

# General Keyword Arguments
 - `finput = default_coupling(submodels...)`: Coupling function for the combined system.
    Has the function signature: `y₁, y₂, ... yₙ = finput(dxᵢ, xᵢ, pᵢ, t)`.  Default coupling
    functions are defined for various combinations of built-in models.
 - `fparam = default_parameter_function(finput)`: Parameter function for the coupled system.
    Has the function signature: `p₁, p₂, ... pₙ, pₐ = fparam(p, t)`

# Keyword Arguments for Automatically Determining Jacobian Sparsity
 - `neval=3`: Total number of evaluation points (for computing jacobian properties)
 - `dx0=[rand(sum(nstate) for ieval = 1:neval]`: Evaluation point rates
 - `x0=[rand(sum(nstate) for ieval = 1:neval]`: Evaluation point states
 - `p0=fill(default_parameters(finput), neval)`: Evaluation point parameters
 - `t0=[rand() for ieval=1:neval]`: Evaluation point times
 - `symbolic=true`: Symbolically find jacobian sparsity?
 - `autodiff=true`: Use automatic differentiation to evaluate jacobians?
 - `fdtype=Val(:forward)`: Finite differencing type to use if `autodiff=false`
 - `ztol=nothing`: Tolerance for determining jacobian sparsity (the default includes structural zeros)
 - `atol=0.0`: Absolute tolerance for determining rate jacobian properties (see `isapprox`)
 - `rtol=0.0`: Relative tolerance for determining rate jacobian properties (see `isapprox`)

 # Keyword Arguments for Providing Jacobian Sparsity (and other properties)
 - `rate_jacobian_sparsity=nothing`: Coupled system rate jacobian sparsity
 - `state_jacobian_sparsity=nothing`: Coupled system state jacobian sparsity
 - `rate_jacobian_colorvec=nothing`: Coupled system rate jacobian color vector
 - `state_jacobian_colorvec=nothing`: Coupled system state jacobian color vector
 - `identity_mass_matrix=nothing`: Indicates whether the mass matrix is the identity matrix
 - `constant_mass_matrix=nothing`: Indicates whether the mass matrix does not change.  Note
    that when `neval=1` this property defaults to `false` unless `identity_mass_matrix==true`
 - `mass_matrix=spzeros(sum(nstate),sum(nstate))`: Storage for the mass matrix.  Note that
    this matrix is only used if `constant_mass_matrix=true`
"""
function CoupledModel(submodels, parameters, nstate=number_of_states.(submodels);
    finput=default_coupling(submodels...), fparam=default_parameter_function(finput), neval=3,
    dx0=[rand(sum(nstate)) for ieval in 1:neval], x0=[rand(sum(nstate)) for ieval in 1:neval],
    p0=fill(parameters, neval), t0=[rand() for ieval in 1:neval],
    symbolic=true, autodiff=true, fdtype=Val(:forward), ztol=nothing, atol=0.0, rtol=0.0,
    rate_jacobian_sparsity=nothing, rate_jacobian_colorvec=nothing,
    state_jacobian_sparsity=nothing, state_jacobian_colorvec=nothing,
    identity_mass_matrix=nothing, constant_mass_matrix=nothing,
    mass_matrix=spzeros(sum(nstate), sum(nstate)))

    # check inputs
    @assert length(submodels) == length(nstate)
    @assert neval == length(dx0) == length(x0) == length(p0) == length(t0)

    # find indices corresponding to each submodel
    nsub = length(submodels)
    lastindex = vcat(0, cumsum(nstate)...)
    indices = [lastindex[isub]+1:lastindex[isub+1] for isub = 1:nsub]

    # set default rate jacobian color vector
    if isnothing(rate_jacobian_colorvec)
        if isnothing(rate_jacobian_sparsity)
            rate_jacobian_colorvec = 1:length(x0[1])
        else
            rate_jacobian_colorvec = SparseDiffTools.matrix_colors(rate_jacobian_sparsity)
        end
    end

    # set default state jacobian color vector
    if isnothing(state_jacobian_colorvec)
        if isnothing(state_jacobian_sparsity)
            state_jacobian_colorvec = 1:length(x0[1])
        else
            state_jacobian_colorvec = SparseDiffTools.matrix_colors(state_jacobian_sparsity)
        end
    end

    # find rate jacobian sparsity and compute sample jacobians
    if symbolic
        # compute jacobian sparsity (symbolically)
        if isnothing(rate_jacobian_sparsity)
            rate_jacobian_sparsity = symbolic_coupled_rate_sparsity_pattern(dx0[1], x0[1], p0[1], t0[1],
                submodels, finput, fparam, indices; ztol=ztol)
        end
        # compute color vector
        if isnothing(rate_jacobian_colorvec)
            rate_jacobian_colorvec = SparseDiffTools.matrix_colors(rate_jacobian_sparsity)
        end
        # sample jacobians (if necessary)
        if isnothing(identity_mass_matrix) || isnothing(constant_mass_matrix)
            rate_jacobs = [zeros(sum(nstate), sum(nstate)) for ieval = 1:neval]
            for ieval = 1:neval
                if autodiff
                    autodiff_coupled_rate_jacobian!(rate_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=rate_jacobian_colorvec, sparsity=rate_jacobian_sparsity)
                else
                    finitediff_coupled_rate_jacobian!(rate_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=rate_jacobian_colorvec, sparsity=rate_jacobian_sparsity,
                        fdtype=fdtype)
                end
            end
        end
    else
        # sample jacobians
        if isnothing(identity_mass_matrix) || isnothing(constant_mass_matrix) || isnothing(rate_jaocbian_sparsity)
            rate_jacobs = [zeros(sum(nstate), sum(nstate)) for ieval = 1:neval]
            for ieval = 1:neval
                if autodiff
                    autodiff_coupled_rate_jacobian!(rate_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=rate_jacobian_colorvec, sparsity=rate_jacobian_sparsity)
                else
                    finitediff_coupled_rate_jacobian!(rate_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=rate_jacobian_colorvec, sparsity=rate_jacobian_sparsity,
                        fdtype=fdtype)
                end
            end
        end
        # compute jacobian sparsity
        if isnothing(rate_jacobian_sparsity)
            rate_jacobian_sparsity = numeric_sparsity_pattern(rate_jacobs...; ztol=ztol)
        end
        # compute color vector
        if isnothing(rate_jacobian_colorvec)
            rate_jacobian_colorvec = SparseDiffTools.matrix_colors(rate_jacobian_sparsity)
        end
    end

    # find state jacobian sparsity
    if symbolic
        # compute jacobian sparsity (symbolically)
        if isnothing(state_jacobian_sparsity)
            state_jacobian_sparsity = symbolic_coupled_state_sparsity_pattern(dx0[1], x0[1], p0[1], t0[1],
                submodels, finput, fparam, indices; ztol=ztol)
        end
        # compute color vector
        if isnothing(state_jacobian_colorvec)
            state_jacobian_colorvec = SparseDiffTools.matrix_colors(state_jacobian_sparsity)
        end
    else
        if isnothing(state_jacobian_sparsity)
            # sample jacobians
            state_jacobs = [zeros(sum(nstate), sum(nstate)) for ieval = 1:neval]
            for ieval = 1:neval
                if autodiff
                    autodiff_coupled_state_jacobian!(state_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=state_jacobian_colorvec, sparsity=state_jacobian_sparsity)
                else
                    finitediff_coupled_state_jacobian!(state_jacobs[ieval], dx0[ieval],
                        x0[ieval], p0[ieval], t0[ieval], submodels, finput, fparam, indices;
                        colorvec=state_jacobian_colorvec, sparsity=state_jacobian_sparsity,
                        fdtype=fdtype)
                end
            end
            # compute jacobian sparsity
            state_jacobian_sparsity = numeric_sparsity_pattern(state_jacobs...; ztol=ztol)
        end
        # compute color vector
        if isnothing(state_jacobian_colorvec)
            state_jacobian_colorvec = SparseDiffTools.matrix_colors(state_jacobian_sparsity)
        end
    end

    # check if the mass matrix is the identity matrix
    if isnothing(identity_mass_matrix)
        identity_mass_matrix = isidentity(rate_jacobs...; atol=atol, rtol=rtol)
    end

    # check if the mass matrix is constant
    if identity_mass_matrix
        constant_mass_matrix = true
    elseif length(rate_jacobs) == 1
        constant_mass_matrix = false
    else
        if isnothing(constant_mass_matrix)
            constant_mass_matrix = isconstant(rate_jacobs...; atol=atol, rtol=rtol)
        end
    end

    # save mass matrix
    mass_matrix .= -rate_jacobs[1]

    return CoupledModel(submodels, finput, fparam, indices, rate_jacobian_sparsity,
        state_jacobian_sparsity, rate_jacobian_colorvec, state_jacobian_colorvec,
        identity_mass_matrix, constant_mass_matrix, mass_matrix)
end

"""
    number_of_states(coupled_model::CoupledModel)

Return the number of state variables in a coupled model
"""
number_of_states(coupled_model::CoupledModel) = coupled_model.indices[end][end]

"""
    state_indices(coupled_model::CoupledModel)

Return indices for accessing submodel states.
"""
state_indices(coupled_model::CoupledModel) = coupled_model.indices

"""
    separate_states(x, coupled_model::CoupledModel)

Return the states corresponding to each submodel.
"""
separate_states(x, coupled_model::CoupledModel) = view.(Ref(x), coupled_model.indices)

"""
    default_parameter_function(finput)

Return the default parameter function for a given coupling function.
"""
default_parameter_function(finput) = (p, t) -> p

"""
    (coupled_model::CoupledModel)(dx, x, p, t)

Compute the residual (out-of-place) for a coupled model.
"""
(coupled_model::CoupledModel)(dx, x, p, t) = residual(coupled_model, dx, x, p, t)

"""
    (coupled_model::CoupledModel)(dx, x, p, t)

Compute the residual (in-place) for a coupled model.
"""
(coupled_model::CoupledModel)(resid, dx, x, p, t) = residual!(resid, coupled_model, dx, x, p, t)

"""
    residual(coupled_model, dx, x, p, t)

Compute the residual for a coupled model.
"""
function residual(coupled_model, dx, x, p, t)

    resid = zeros(eltype(dx), length(dx))

    return residual!(resid, coupled_model, dx, x, p, t)
end

"""
    residual!(resid, coupled_model, dx, x, p, t)

Compute the residual for a coupled model.
"""
function residual!(resid, coupled_model, dx, x, p, t)

    @unpack fresid, finput, fparam, indices = coupled_model

    return coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)
end

"""
    rate_jacobian(coupled_model, dx, x, p, t; kwargs...)

Compute the jacobian of the residual with respect to the state rates.
"""
function rate_jacobian(coupled_model, dx, x, p, t; kwargs...)

    jacob = spzeros(eltype(dx), length(dx), length(dx))

    return rate_jacobian!(jacob, coupled_model, dx, x, p, t; kwargs...)
end

"""
    rate_jacobian!(jacob, coupled_model, dx, x, p, t; kwargs...)

Compute the jacobian of the residual with respect to the state rates.
"""
function rate_jacobian!(jacob, coupled_model, dx, x, p, t; autodiff=true, fdtype=Val(:forward))

    @unpack fresid, finput, fparam, indices = coupled_model

    colorvec = coupled_model.rate_jacobian_colorvec
    sparsity = coupled_model.rate_jacobian_sparsity

    if autodiff
        autodiff_coupled_rate_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec, sparsity)
    else
        finitediff_coupled_rate_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec, sparsity, fdtype)
    end

    return jacob
end

"""
    state_jacobian(coupled_model, dx, x, p, t; kwargs...)

Compute the jacobian of the residual with respect to the states.
"""
function state_jacobian(coupled_model, dx, x, p, t; kwargs...)

    jacob = spzeros(eltype(dx), length(dx), length(dx))

    return state_jacobian!(jacob, coupled_model, dx, x, p, t; kwargs...)
end

"""
    state_jacobian!(jacob, coupled_model, dx, x, p, t; kwargs...)

Compute the jacobian of the residual with respect to the state rates.
"""
function state_jacobian!(jacob, coupled_model, dx, x, p, t; autodiff=true, fdtype=Val(:forward))

    @unpack fresid, finput, fparam, indices = coupled_model

    colorvec = coupled_model.state_jacobian_colorvec
    sparsity = coupled_model.state_jacobian_sparsity

    if autodiff
        autodiff_coupled_state_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec, sparsity)
    else
        finitediff_coupled_state_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec, sparsity, fdtype)
    end

    return jacob
end

function coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    # separate residuals
    residuals = view.(Ref(resid), indices)

    # separate rates
    rates = view.(Ref(dx), indices)

    # separate states
    states = view.(Ref(x), indices)

    # update parameters
    parameters = fparam(p, t)

    # compute inputs
    inputs = finput(rates, states, parameters, t)

    # compute residuals (in-place)
    for i in eachindex(fresid)
        fresid[i](residuals[i], rates[i], states[i], inputs[i], parameters[i], t)
    end

    # return modified residual vector
    return resid
end

function finitediff_coupled_rate_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec=1:length(dx), sparsity=nothing, fdtype=Val(:forward))

    f! = (resid, dx) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    jac_cache = FiniteDiff.JacobianCache(dx, fdtype; colorvec, sparsity)

    return FiniteDiff.finite_difference_jacobian!(jacob, f!, dx, jac_cache)
end

function autodiff_coupled_rate_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec=1:length(dx), sparsity=nothing)

    f! = (resid, dx) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    jac_cache = ForwardColorJacCache(f!, dx; colorvec, sparsity)

    return SparseDiffTools.forwarddiff_color_jacobian!(jacob, f!, dx, jac_cache)
end

function symbolic_coupled_rate_sparsity_pattern(dx, x, p, t, fresid, finput, fparam, indices; ztol=ztol)

    f! = (resid, dx) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    return symbolic_sparsity_pattern(f!, similar(dx), dx; ztol=ztol)
end

function finitediff_coupled_state_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec=1:length(dx), sparsity=nothing, fdtype=Val(:forward))

    f! = (resid, x) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    jac_cache = FiniteDiff.JacobianCache(x, fdtype; colorvec, sparsity)

    return FiniteDiff.finite_difference_jacobian!(jacob, f!, x, jac_cache)
end

function autodiff_coupled_state_jacobian!(jacob, dx, x, p, t, fresid, finput, fparam, indices; colorvec=1:length(dx), sparsity=nothing)

    f! = (resid, x) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    jac_cache = ForwardColorJacCache(f!, x; colorvec, sparsity)

    return SparseDiffTools.forwarddiff_color_jacobian!(jacob, f!, x, jac_cache)
end

function symbolic_coupled_state_sparsity_pattern(dx, x, p, t, fresid, finput, fparam, indices; ztol=ztol)

    f! = (resid, x) -> coupled_residual!(resid, dx, x, p, t, fresid, finput, fparam, indices)

    return symbolic_sparsity_pattern(f!, similar(x), x; ztol=ztol)
end

function symbolic_sparsity_pattern(f!, output, inputs...; ztol=nothing)
    # find sparsity pattern
    sparsity_pattern = spzeros(length(output), length(inputs[1]))
    for input in inputs
        J = Symbolics.jacobian_sparsity(f!,output,input)
        sparsity_pattern .+= abs.(J)
    end
    sparsity_pattern = sparse(sparsity_pattern)
    # apply tolerance (if desired)
    isnothing(ztol) || droptol!(sparsity_pattern, ztol)
    # return sparsity pattern
    return sparsity_pattern
end

function numeric_sparsity_pattern(jacobs...; ztol=nothing)
    # find sparsity pattern
    sparsity_pattern = similar(jacobs[1]) .= 0
    for jacob in jacobs
        sparsity_pattern .+= abs.(jacob)
    end
    sparsity_pattern = sparse(sparsity_pattern)
    # apply tolerance (if desired)
    isnothing(ztol) || droptol!(sparsity_pattern, ztol)
    # return sparsity pattern
    return sparsity_pattern
end

function isidentity(jacobs...; atol=0, rtol=0)
    result = true
    for jacob in jacobs
        result = result && isapprox(jacob, I)
    end
    return result
end

function isconstant(jacobs...; atol=0, rtol=0)
    result = true
    for jacob in jacobs
        result = result && isapprox(jacobs[1], jacob)
    end
    return result
end