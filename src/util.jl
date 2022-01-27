# --- extract values from Val() objects --- #

val(x) = x
val(::Val{x}) where x = x

# --- automatic differentiation --- #

# out-of-place
function autodiff_jacobian_func(f, ivar)
    (args...) -> begin
        x = args[ivar]
        fy = (x) -> f(Base.setindex(args, x, ivar)...)
        ForwardDiff.jacobian(fy, x)
    end
end

# in-place
function autodiff_jacobian_func(f, ny, ivar)
    (J, args...) -> begin
        x = args[ivar]
        y = similar(J, val(ny))
        fy = (y, x) -> f(y, Base.setindex(args, x, ivar)...)
        ForwardDiff.jacobian!(J, fy, y, x)
    end
end

# out-of-place
function autodiff_derivative_func(f, ivar)
    (args...) -> begin
        x = args[ivar]
        fy = (x) -> f(Base.setindex(args, x, ivar)...)
        ForwardDiff.derivative(fy, x)
    end
end

# in-place
function autodiff_derivative_func(f, ny, ivar)
    (dy, args...) -> begin
        x = args[ivar]
        y = similar(dy, val(ny))
        fy = (y, x) -> f(y, Base.setindex(args, x, ivar)...)
        ForwardDiff.derivative!(dy, fy, y, x)
    end
end

# --- block diagonal --- #

function block_diagonal(Jii)
    
    J = zeros(sum(size.(Jii, 1)), sum(size.(Jii, 2)))

    i = j = 1
    for iblock = 1:N
        nx, ny = size(Jii[iblock])
        ix, iy = i:i+nx-1, j:j+ny-1
        J[ix, iy] = Jii[iblock]
        i += nx; j += ny
    end
    
    return J
end

block_diagonal(Jii::NTuple{N,StaticMatrix}) where N = static_block_diagonal(Jii)

@generated function static_block_diagonal(Jii::NTuple{N,StaticMatrix}) where N

    # assemble first row of the jacobian matrix
    expr = quote
        J = Jii[1]
    end
    for j = 2:N
        if 1 == j
            expr = quote
                $expr
                Jij = Jii[1]
            end
        else
            expr = quote
                $expr
                nx = Size(Jii[1])[1]
                ny = Size(Jii[$j])[2]
                sz = Size(nx, ny)
                Jij = zero(similar_type(Jii[1], sz))
            end
        end
        expr = quote
            $expr
            J = hcat(J, Jij)
        end
    end

    # add remaining rows of the jacobian matrix
    for i = 2:N
        if i == 1
            expr = quote
                $expr
                Ji = Jii[$i]
            end
        else
            expr = quote
                $expr
                nx = Size(Jii[$i])[1]
                ny = Size(Jii[1])[2]
                sz = Size(nx, ny)
                Ji = zero(similar_type(Jii[$i], sz))
            end
        end
        for j = 2:N
            if i == j
                expr = quote
                    $expr
                    Jij = Jii[$i]
                end
            else
                expr = quote
                    $expr
                    nx = Size(Jii[$i])[1]
                    ny = Size(Jii[$j])[2]
                    sz = Size(nx, ny)
                    Jij = zero(similar_type(Jii[$i], sz))
                end
            end
            expr = quote
                $expr
                Ji = hcat(Ji, Jij)
            end
        end
        expr = quote
            $expr
            J = vcat(J, Ji)
        end
    end

    return expr
end

# --- eigenmode correlation --- #

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

