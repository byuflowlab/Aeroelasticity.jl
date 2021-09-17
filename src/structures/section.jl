"""
    TypicalSection <: UnsteadyModel

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct TypicalSection <: UnsteadyModel end

"""
    TypicalSection()

Initialize an object of type [`TypicalSection`](@ref)
"""
TypicalSection()

# --- Traits --- #

number_of_states(::Type{TypicalSection}) = 4
number_of_inputs(::Type{TypicalSection}) = 2
number_of_parameters(::Type{TypicalSection}) = 5

inplaceness(::Type{TypicalSection}) = OutOfPlace()

rate_jacobian_type(::Type{TypicalSection}) = Constant()
state_jacobian_type(::Type{TypicalSection}) = Constant()
input_jacobian_type(::Type{TypicalSection}) = Invariant()
parameter_jacobian_type(::Type{TypicalSection}) = Linear()
time_gradient_type(::Type{TypicalSection}) = Zeros()

# --- Methods --- #

function get_residual(::TypicalSection, dx, x, y, p, t)
    # extract rates
    dh, dθ, dhdot, dθdot = dx
    # extract states
    h, θ, hdot, θdot = x
    # extract inputs
    L, M = y
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate state rates
    return section_residual(dh, dθ, dhdot, dθdot, h, θ, hdot, θdot, L, M,
        kh, kθ, m, Sθ, Iθ)
end

# --- Performance Overloads --- #

function get_rate_jacobian(::TypicalSection, p)
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return section_rate_jacobian(m, Sθ, Iθ)
end

function get_state_jacobian(::TypicalSection, p)
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return section_state_jacobian(kh, kθ)
end

function get_input_jacobian(::TypicalSection)
    # return jacobian
    return section_input_jacobian()
end

function get_parameter_jacobian(::TypicalSection, dx, x, y, p, t)
    # extract rates
    dh, dθ, dhdot, dθdot = dx
    # extract states
    h, θ, hdot, θdot = x
    # extract inputs
    L, M = y
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return section_parameter_jacobian(h, θ, dhdot, dθdot)
end

# --- Convenience Methods --- #

function set_states!(x, model::TypicalSection; h, theta, hdot, thetadot)

    x[1] = h
    x[2] = theta
    x[3] = hdot
    x[4] = thetadot

    return x
end

function set_inputs!(y, model::TypicalSection; L, M)

    y[1] = L
    y[2] = M

    return y
end

function set_parameters!(p, model::TypicalSection; kh, ktheta, m, Stheta, Itheta)

    p[1] = kh
    p[2] = ktheta
    p[3] = m
    p[4] = Stheta
    p[5] = Itheta

    return p
end

function separate_states(model::TypicalSection, x)

    return (h = x[1], theta = x[2], hdot = x[3], thetadot = x[4])
end

function separate_inputs(model::TypicalSection, y)

    return (L = y[1], M = y[2])
end

function separate_parameters(model::TypicalSection, p)

    return (kh = p[1], ktheta = p[2], m = p[3], Stheta = p[4], Itheta = p[5])
end

# --- Plotting --- #

@recipe function f(model::TypicalSection, dx, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xticks --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    h = x[1]
    θ = x[2]

    xplot = [-0.5*cos(θ),  0.5*cos(θ)]
    yplot = [ 0.5*sin(θ)-h, -0.5*sin(θ)-h]

    xplot, yplot
end

# --- Internal Methods for this Model --- #

function section_residual(dh, dθ, dhdot, dθdot, h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)
    r1 = dh - hdot
    r2 = dθ - θdot
    r3 = m*dhdot + Sθ*dθdot + kh*h + L
    r4 = Sθ*dhdot + Iθ*dθdot + kθ*θ - M
    return SVector(r1, r2, r3, r4)
end

section_rate_jacobian(m, Sθ, Iθ) = @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]

section_state_jacobian(kh, kθ) = @SMatrix [0 0 -1 0; 0 0 0 -1; kh 0 0 0; 0 kθ 0 0]

section_input_jacobian() = @SMatrix [0 0; 0 0; 1 0; 0 -1]

function section_parameter_jacobian(h, θ, dhdot, dθdot)
    return @SMatrix [0 0 0 0 0; 0 0 0 0 0; h 0 dhdot dθdot 0; 0 θ 0 dhdot dθdot]
end

# --- Internal Methods for Couplings --- #

# section velocities for coupling with 2D aerodynamic models
function section_velocities(U, θ, hdot, θdot)

    u = U
    v = U*θ + hdot
    ω = θdot

    return SVector(u, v, ω)
end

# section accelerations for coupling with 2D aerodynamic models
function section_accelerations(dhdot, dθdot)
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    return SVector(udot, vdot, ωdot)
end

# derivatives
section_velocities_U(θ) = SVector(1, θ, 0)
section_velocities_θ(U) = SVector(0, U, 0)
section_velocities_hdot() = SVector(0, 1, 0)
section_velocities_θdot() = SVector(0, 0, 1)
section_accelerations_dhdot() = SVector(0, 1, 0)
section_accelerations_dθdot() = SVector(0, 0, 1)
