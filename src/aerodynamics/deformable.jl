"""
    PetersDeformable{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``u, \\omega,
\\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0``
"""
struct PetersDeformable{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #

"""
    PetersDeformable{N,TF=Float64}()

Initialize an object of type `PetersDeformable` which has `N` aerodynamic degrees of
freedom.
"""
PetersDeformable{N}() where N = PetersDeformable{N,Float64}()

function PetersDeformable{N,TF}() where {N,TF}

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

    return PetersDeformable(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #

number_of_states(::Type{PetersDeformable{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:PetersDeformable}) = 4
number_of_parameters(::Type{<:PetersDeformable}) = 4

inplaceness(::Type{<:PetersDeformable}) = OutOfPlace()

rate_jacobian_type(::Type{<:PetersDeformable}) = Invariant()
state_jacobian_type(::Type{<:PetersDeformable}) = Linear()
input_jacobian_type(::Type{<:PetersDeformable}) = Nonlinear()
parameter_jacobian_type(::Type{<:PetersDeformable}) = Nonlinear()
time_gradient_type(::Type{<:PetersDeformable}) = Zeros()

# --- Methods --- #

function get_residual(model::PetersDeformable{N,TF,SV,SA}, dx, x, y, p, t) where {N,TF,SV,SA}
    # extract rates
    dλ = SVector{N}(dx)
    # extract states
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    Abar = model.A
    cbar = model.c
    # calculate rates
    return peters_residual(dλ, λ, u, ω, vdot, ωdot, a, b, Abar, cbar)
end

# --- Performance Overloads --- #

# TODO

# --- Convenience Functions --- #

# TODO

# --- Internal Methods for Model --- #

function peters_residual(dλ, λ, u, ω, vdot, ωdot, a, b, Abar, cbar)

    return Abar*dλ - cbar*(vdot + u*ω + (b/2-a*b)*ωdot) + u/b*λ
end

# --- Internal Methods for Couplings --- #


function peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    # induced flow velocity
    λ0 = 1/2 * bbar'*λ

    # no reversed flow
    f = 1

    # adjust section velocities and accelerations
    u = u # chordwise velocity
    vt = v + (b/2 - a*b)*ω - λ0 - u*α0  # normal velocity at 3/4 chord
    ω = ω # angular velocity
    udot = udot # chordwise acceleration
    vdot = vdot - u*ω - a*b*ωdot # normal acceleration at midchord
    ωdot = ωdot # angular acceleration

    # compute normal force
    N = a0*f*ρ*b*u*vt + 2*pi*ρ*b^2*u*ω + pi*ρ*b^2*vdot

    # compute axial force
    A = -a0*f*ρ*b*vt^2

    # moment at the reference location
    M = a0/2*f*ρ*b^2*u*vt - pi*ρ*b^4/8*ωdot - pi*ρ*b^3/2*u*ω + a*b*N

    # return loads at the reference location
    return SVector(N, A, M)
end


function peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    # induced flow velocity
    λ0 = 1/2 * bbar'*λ

    # no reversed flow
    f = 1

    # adjust section velocities and accelerations
    u = u # chordwise velocity
    vt = v + (b/2 - a*b)*ω - λ0 - u*α0  # normal velocity at 3/4 chord
    ω = ω # angular velocity
    udot = udot # chordwise acceleration
    vdot = vdot - u*ω - a*b*ωdot # normal acceleration at midchord
    ωdot = ωdot # angular acceleration

    # compute generalized loads
    tmpsum = u0/b*(f*sum(k -> k*h[k], 1:3:Nh) + sum(k -> k*h[k], 2:4:Nh))
    Y = f(v0 + h0dot - λ) + tmpsum

    # compute normal force
    -2*pi*ρ*b^2*N = - u0/b*(vt + h1dot) -1/2*(-vdot + h0ddot - h2ddot/2) - u0/b*Y  - 1/2*(u0dot/b*h1 - u0/b*vt)

    # compute axial force
    tmpsum = sum(k -> h[k]*((hddot[k-1] - hddot[k+1])/4 + u0/b*k*hdot[k] + 1/2*u0dot/b*k*h[k]))
    2*pi*ρ*b*A = -(v0 + h0dot - λ0)*Y + tmpsum + (-vdot/2 + h0ddot/4 + u0/b*vt/2)*h1 + 1/4*v1dot*h2

    # compute moment at reference point
    2*pi*ρ*b*M = -1/16*(v1dot + h1ddot - h3ddot) + 1/2*u0/b*(u0 + h0dot - λ0) -
        1/2*u0/b*h2dot + 1/2*u0^2/b^2*h1 - 1/4*u0dot/b*h2 - a*b*N

    # return loads
    return SVector(N, A, M)
end
