"""
    Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}}

Two-dimensional aerodynamic model based on Peters' finite state model with `N` state
variables, inputs ``u, \\omega, \\dot{v}, \\dot{\\omega}`` and parameters
``a, b, a_0, \\alpha_0, c_{d,0}, c_{m,0}``
"""
struct Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}}
    A::TA
    b::TV
    c::TV
end

"""
    Peters{N, TF=Float64}()

Initialize a model with type [`Peters`](@ref) with `N` state variables
"""
Peters()

Peters{N}() where N = Peters{N,Float64}()

function Peters{N,TF}() where {N,TF}

    A, b, c = peters_constants(N, TF)

    return Peters(A, b, c)
end

# residual function
function (peters::Peters)(resid, dx, x, y, p, t)

    # extract constants
    Abar = peters.A
    cbar = peters.c

    # create static versions of dλ and λ
    dλ = similar_type(cbar, eltype(dx))(dx)
    λ = similar_type(cbar, eltype(x))(x)

    # extract inputs
    u, ω, vdot, ωdot = y

    # extract parameters
    a, b, a0, α0, cd0, cm0 = p

    # compute and return residual
    resid .= Abar*dλ - cbar*(vdot + u*ω + (b/2-a*b)*ωdot) + u/b*λ
end

number_of_states(::Peters{N,TF,TV,TA}) where {N,TF,TV,TA} = N

number_of_parameters(::Peters) = 6

# --- Internal Methods --- #

# function for defining model constants
function peters_constants(N, TF)

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return SMatrix{N,N}(A), SVector{N}(b), SVector{N}(c)
end

# aerodynamic loads per unit span
function peters_loads(a, b, a0, α0, cd0, cm0, rho, beta, bbar, u, v, ω, vdot, ωdot, λ)
    # circulatory load factor
    tmp1 = a0*rho*u*b
    # non-circulatory load factor
    tmp2 = pi*rho*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # normal force at the reference point
    N = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at the reference point
    A = -a0*rho*b*(v + d*ω - λ0 - u*α0)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + 2*rho*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction (use constant beta to preserve linearity)
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag (note that Prandtl-Glauert correction applies only to pressure)
    A += rho*b*u^2*cd0
    # return loads
    return SVector(N, A, M)
end