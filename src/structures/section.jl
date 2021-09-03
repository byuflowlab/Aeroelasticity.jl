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
rate_jacobian_type(::Type{TypicalSection}) = Identity()
state_jacobian_type(::Type{TypicalSection}) = Constant()
input_jacobian_type(::Type{TypicalSection}) = Constant()
parameter_jacobian_type(::Type{TypicalSection}) = Linear()

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

function get_state_jacobian(::TypicalSection, p)
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return -section_state_jacobian(kh, kθ, m, Sθ, Iθ)
end

function get_input_jacobian(::TypicalSection, p)
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return -section_input_jacobian(kh, kθ, m, Sθ, Iθ)
end

function get_parameter_jacobian(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return -section_parameter_jacobian(h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)
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

@recipe function f(model::TypicalSection, x, y, p, t)

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

# --- Internal Methods --- #

# state rate residual
function section_residual(dh, dθ, dhdot, dθdot, h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)
    dx = SVector(dh, dθ, dhdot, dθdot)
    xdot = section_rates(h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)
    return dx - xdot
end

# state rates
function section_rates(h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)
    dh = hdot
    dθ = θdot
    dhdot = (-Iθ*kh*h + Sθ*kθ*θ + Iθ*L + Sθ*M) / (m*Iθ - Sθ^2)
    dθdot = ( Sθ*kh*h - m*kθ*θ + Sθ*L + m*M) / (m*Iθ - Sθ^2)
    return SVector(dh, dθ, dhdot, dθdot)
end

# jacobian of state rates wrt states
function section_state_jacobian(kh, kθ, m, Sθ, Iθ)
    tmp = m*Iθ - Sθ^2
    @SMatrix [0 0 1 0; 0 0 0 1; -Iθ*kh/tmp Sθ*kθ/tmp 0 0; Sθ*kh/tmp -m*kθ/tmp 0 0]
end

# jacobian of state rates wrt inputs
function section_input_jacobian(kh, kθ, m, Sθ, Iθ)
    tmp = m*Iθ - Sθ^2
    return @SMatrix [0 0; 0 0; Iθ/tmp Sθ/tmp; Sθ/tmp m/tmp]
end

# jacobian of state rates wrt parameters
function section_parameter_jacobian(h, θ, hdot, θdot, L, M, kh, kθ, m, Sθ, Iθ)

    tmp = 1/(m*Iθ - Sθ^2)
    tmp_m = -Iθ/(m*Iθ - Sθ^2)^2
    tmp_S = 2*Sθ/(m*Iθ - Sθ^2)^2
    tmp_Iθ = -m/(m*Iθ - Sθ^2)^2

    dhdot = (-Iθ*kh*h + Sθ*kθ*θ + Iθ*L + Sθ*M) * tmp
    dhdot_kh = -Iθ*h*tmp
    dhdot_kθ = Sθ*θ*tmp
    dhdot_m = dhdot*tmp_m
    dhdot_S = (kθ*θ + M) * tmp + dhdot*tmp_S
    dhdot_Iθ = (-kh*h + L) * tmp + dhdot*tmp_Iθ

    dθdot = ( Sθ*kh*h - m*kθ*θ + Sθ*L + m*M) * tmp
    dθdot_kh = Sθ*h*tmp
    dθdot_kθ = -m*θ*tmp
    dθdot_m = (-kθ*θ + M)*tmp + dθdot*tmp_m
    dθdot_S = (kh*h + L) * tmp + dθdot*tmp_S
    dθdot_Iθ = dθdot*tmp_Iθ

    return @SMatrix [
        0 0 0 0 0;
        0 0 0 0 0;
        dhdot_kh dhdot_kθ dhdot_m dhdot_S dhdot_Iθ;
        dθdot_kh dθdot_kθ dθdot_m dθdot_S dθdot_Iθ;
        ]
end

function section_velocities(U, θ, hdot, θdot)

    u = U
    v = U*θ + hdot
    ω = θdot

    return SVector(u, v, ω)
end

section_velocities_U(U, θ, hdot, θdot) = SVector(1, θ, 0)
section_velocities_θ(U, θ, hdot, θdot) = SVector(0, U, 0)
section_velocities_hdot(U, θ, hdot, θdot) = SVector(0, 1, 0)
section_velocities_θdot() = SVector(0, 0, 1)

function section_accelerations(dhdot, dθdot)
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    return SVector(udot, vdot, ωdot)
end

section_accelerations_dhdot() = SVector(0, 1, 0)
section_accelerations_dθdot() = SVector(0, 0, 1)
