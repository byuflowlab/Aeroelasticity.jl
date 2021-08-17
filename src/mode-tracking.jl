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
