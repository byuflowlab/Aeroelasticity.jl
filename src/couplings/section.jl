# --- TypicalSection + PetersFiniteState --- #
inplace_input(::PetersFiniteState, ::TypicalSection) = false
has_input_mass_matrix(::PetersFiniteState, ::TypicalSection) = true
constant_input_mass_matrix(::PetersFiniteState, ::TypicalSection) = false
defined_input_state_jacobian(::PetersFiniteState, ::TypicalSection) = true

function get_input_mass_matrix(aero::PetersFiniteState, stru::TypicalSection, u, p, t)
    # extract model constants
    cbar = model.c
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # create zero row vector with length nλ
    zλ = zero(aero.c')
    # construct submatrices
    Mds = @SMatrix [0 0 0 0; 0 0 -1 0; 0 0 0 -1]
    Mda = @SMatrix vcat(zλ, zλ, zλ, zλ)
    Mrs = hcat(
        @SVector zeros(2),
        @SVector zeros(2),
        peters_loads_hddot(b, ρ)
        peters_loads_θddot(a, b, ρ)
        )
    Mra = peters_loads_λdot(zλ)
    # assemble mass matrix
    return vcat(hcat(Mds, Mda), hcat(Mrs, Mra))
end

function get_inputs(aero::PetersFiniteState{N,TF}, stru::TypicalSection,
    u, p, t) where {N,TF}
    # indices for extracting state variables
    iq = SVector{4}(1:4)
    iλ = SVector{N}(5:4+N)
    # separate aerodynamic and structural states
    q = u[iq]
    λ = u[iλ]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # calculate aerodynamic loads
    L, M = peters_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)
    # return portion of inputs that is not dependent on the state rates
    return SVector(θdot, 0, 0, L, M)
end

function get_input_state_jacobian(aero::PetersFiniteState, stru::TypicalSection,
    u, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # create zero row vector with length nλ
    zλ = zero(aero.bbar')
    # compute jacobian sub-matrices
    Jds = @SMatrix [0 0 0 1; 0 0 0 0; 0 0 0 0]
    Jda = vcat(zλ, zλ, zλ)
    Jrs = hcat(
        peters_loads_h(b, U, ρ),
        peters_loads_θ(b, U, ρ),
        peters_loads_hdot()
        peters_loads_θdot(a, b, U, ρ)
        )
    Jra = peters_loads_λ(b, ρ, bbar)
    # return jacobian
    return vcat(hcat(Jds, Jda), hcat(Jrs, Jra))
end

# TODO: Parameter jacobian
