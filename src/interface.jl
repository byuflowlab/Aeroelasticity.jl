"""
    linearize(model, x, p; autodiff=true, fdtype=Val(:forward))
    linearize(model, dx, x, p, t; autodiff=true, fdtype=Val(:forward))

Calculate the jacobian of the residual function for coupled model `model` with respect to
the state variables and their rates.
"""
linearize

linearize(model::CoupledModel, x, p; kwargs...) = linearize(model, FillArrays.Zeros(x), x, p, 0; kwargs...)

function linearize(model::CoupledModel, dx, x, p, t; autodiff=true, fdtype=Val(:forward))

    K = state_jacobian(model, dx, x, p, t; autodiff=autodiff, fdtype=fdtype)
    M = rate_jacobian(model, dx, x, p, t; autodiff=autodiff, fdtype=fdtype)

    return K, M
end

"""
    dense_eigen(K, M; kwargs...)

Return the eigenvalues, left eigenvector matrix, and right eigenvector matrix
corresponding to the model.
"""
function dense_eigen(K, M)
    A = Array(K)
    B = Array(-M)
    E = eigen(A, B) # eigenvalue decomposition
    λ = eigvals(E) # eigenvalues
    V = eigvecs(E) # right eigenvector matrix
    U = I/(M*V) # left eigenvector matrix
    return λ, U, V
end

"""
    sparse_eigen(K, M; nev=min(20, size(K,1)))

Return the eigenvalues, left eigenvector matrix, and right eigenvector matrix
corresponding to the model. `nev` is the number of eigenvalues to compute
"""
function sparse_eigen(K, M; nev=min(20, size(K,1)))

    # construct linear map
    T = promote_type(eltype(K), eltype(M))
    nx = size(K, 1)
    Kfact = lu(K)
    f! = (b, x) -> ldiv!(b, Kfact, M * x)
    fc! = (b, x) -> mul!(b, M', Kfact' \ x)
    L = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # # compute eigenvalues and eigenvectors
    λ, V, _ = Arpack.eigs(L; nev=min(nx,nev), which=:LM)
    # λ, V = partialeigen(partialschur(L; nev=min(nx, nev), which=LM(), tol=1e-6)[1])

    # sort eigenvalues by magnitude
    perm = sortperm(λ, by=(λ) -> (abs(λ), imag(λ)), rev=true)
    λ = λ[perm]
    V = V[:,perm]

    # eigenvalues are actually -1/λ, no modification necessary for eigenvectors
    λ .= -1 ./ λ

    # also compute left eigenvectors
    U = left_eigenvectors(K, M, λ, V)

    # U = similar(V', complex(T))
    # for iλ in eachindex(λ)
    #     # set shift value
    #     shift = λ[iλ] + 1e-6im
    #     # compute corresponding left eigenvector
    #     F = lu(transpose(K) + shift * transpose(M))
    #     Fmap = LinearMap{complex(T)}((y, x) -> ldiv!(y, F, x), nx, ismutating=true)
    #     _, u = IterativeSolvers.powm(Fmap; inverse = true, shift = shift)
    #     # apply normalization
    #     u ./= (transpose(u)*M*view(V,:,iλ))
    #     # store as a row vector
    #     U[iλ,:] .= u
    # end

    return λ, U, V
end

"""
    correlate_eigenmodes(U, M, V; tracked_modes=1:size(U,1), rtol=1e-3, atol=1e-3)

Correlate modes from a current iteration with those from a previous iteration. Returns a
permutation that may be used to reorder the new modes and the associated corruption indices.
"""
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
    ODEFunction(model::Aeroelasticity.CoupledModel, p0)

Construct an `ODEFunction` for a coupled model.  Note that by using this function you assume
that the governing equations may be expressed in the form `M(x, p, t)*dx = f(x, p, t)`
"""
function SciMLBase.ODEFunction(model::CoupledModel, p0)

    dx = zeros(number_of_states(model))

    f = (resid, x, p, t) -> begin
        residual!(resid, model, dx, x, p, t)
        resid .*= -1
    end

    if model.constant_mass_matrix
        mass_matrix = rate_jacobian!(model.mass_matrix, model, dx, dx, p0, 0.0)
    else
        update_func = (jacob, x, p, t) -> rate_jacobian!(jacob, model, dx, x, p, t)
        mass_matrix = DiffEqArrayOperator(collect(model.mass_matrix); update_func)
    end

    jac = (jacob, x, p, t) -> begin
        state_jacobian!(jacob, model, dx, x, p, t)
        jacob .*= -1
    end

    jac_prototype = model.state_jacobian_sparsity

    sparsity = model.state_jacobian_sparsity

    colorvec = model.state_jacobian_colorvec

    if model.constant_mass_matrix
        odefunc = ODEFunction(f; mass_matrix, jac, jac_prototype, sparsity, colorvec)
    else
        odefunc = ODEFunction(f; mass_matrix, jac)
    end

    return odefunc
end

"""
    DAEFunction(model::Aeroelasticity.CoupledModel)

Construct an `DAEFunction` for a coupled model.
"""
function SciMLBase.DAEFunction(model::CoupledModel)

    f = (resid, dx, x, p, t) -> residual!(resid, model, dx, x, p, t)

    mass_matrix = model.mass_matrix

    jac = (jacob, dx, x, p, gamma, t) -> begin
        # compute jacobian
        state_jacobian!(jacob, model, dx, x, p, t)
        # compute mass matrix
        rate_jacobian!(mass_matrix, model, dx, x, p, t)
        # add mass matrix
        jacob .+= gamma * mass_matrix
    end

    jac_prototype = model.rate_jacobian_sparsity .+ model.state_jacobian_sparsity

    sparsity = model.rate_jacobian_sparsity .+ model.state_jacobian_sparsity

    colorvec = SparseDiffTools.matrix_colors(sparsity)

    return DAEFunction(f)
end
