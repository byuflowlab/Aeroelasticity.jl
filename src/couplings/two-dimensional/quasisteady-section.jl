"""
    couple_models(aero::Steady, stru::TypicalSection)

Create an aerostructural model using a steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` and air density ``\\rho`` as additional parameters.
"""
couple_models(aero::Steady, stru::TypicalSection) = (aero, stru)

"""
    couple_models(aero::QuasiSteady, stru::TypicalSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` and air density ``\\rho`` as additional parameters.
"""
couple_models(aero::QuasiSteady, stru::TypicalSection) = (aero, stru)

# --- Traits --- #

number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Constant()
coupling_parameter_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()

number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Constant()
coupling_parameter_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()

number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Constant()
coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Constant()
coupling_parameter_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::QuasiSteady{0}, stru::TypicalSection, dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v = section_steady_velocities(U, θ)
    # calculate aerodynamic loads
    N, A, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(L, M)
end

function get_coupling_inputs(aero::QuasiSteady{1}, stru::TypicalSection, dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    # calculate aerodynamic loads
    N, A, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(L, M)
end

function get_coupling_inputs(aero::QuasiSteady{2}, stru::TypicalSection, dx, x, p, t)
    # extract state variables
    dh, dθ, dhdot, dθdot = dx
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    # local freestream accelerations
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate aerodynamic loads
    N, A, M = quasisteady2_loads(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(L, M)
end

# --- Performance Overloads --- #

function get_coupling_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, p)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady0_section_state_jacobian(a, b, ρ, a0, U)
end

function get_coupling_parameter_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, dx, x, p, t)

    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v = section_steady_velocities(U, θ)
    u_U, v_U = section_steady_velocities_U(θ)
    # calculate loads
    N_a, A_a, M_a = quasisteady0_loads_a(a, b, ρ, a0, α0, u, v)
    N_b, A_b, M_b = quasisteady0_loads_b(a, b, ρ, a0, α0, u, v)
    N_a0, A_a0, M_a0 = quasisteady0_loads_a0(a, b, ρ, α0, u, v)
    N_α0, A_α0, M_α0 = quasisteady0_loads_α0(a, b, ρ, a0, α0, u, v)

    N_kh, A_kh, M_kh = 0, 0, 0
    N_kθ, A_kθ, M_kθ = 0, 0, 0
    N_m, A_m, M_m = 0, 0, 0
    N_Sθ, A_Sθ, M_Sθ = 0, 0, 0
    N_Iθ, A_Iθ, M_Iθ = 0, 0, 0

    N_u, A_u, M_u = quasisteady0_loads_u(a, b, ρ, a0, α0, u, v)
    N_v, A_v, M_v = quasisteady0_loads_v(a, b, ρ, a0, α0, u, v)
    N_U = N_u * u_U + N_v * v_U
    M_U = M_u * u_U + M_v * v_U

    N_ρ, A_ρ, M_ρ = quasisteady0_loads_ρ(a, b, ρ, a0, α0, u, v)

    # return jacobian
    return @SMatrix [
        N_a N_b N_a0 N_α0 N_kh N_kθ N_m N_Sθ N_Iθ N_U N_ρ;
        M_a M_b M_a0 M_α0 M_kh M_kθ M_m M_Sθ M_Iθ M_U M_ρ
        ]
end

function get_coupling_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, p)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady1_section_state_jacobian(a, b, ρ, a0, U)
end


function get_coupling_parameter_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, dx, x, p, t)

    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    u_U, v_U, ω_U = section_velocities_U(θ)

    # calculate loads
    N_a, A_a, M_a = quasisteady1_loads_a(a, b, ρ, a0, α0, u, v, ω)
    N_b, A_b, M_b = quasisteady1_loads_b(a, b, ρ, a0, α0, u, v, ω)
    N_a0, A_a0, M_a0 = quasisteady1_loads_a0(a, b, ρ, a0, α0, u, v, ω)
    N_α0, A_α0, M_α0 = quasisteady1_loads_α0(a, b, ρ, a0, α0, u, v, ω)

    N_kh, A_kh, M_kh = 0, 0, 0
    N_kθ, A_kθ, M_kθ = 0, 0, 0
    N_m, A_m, M_m = 0, 0, 0
    N_Sθ, A_Sθ, M_Sθ = 0, 0, 0
    N_Iθ, A_Iθ, M_Iθ = 0, 0, 0

    N_u, A_u, M_u = quasisteady1_loads_u(a, b, ρ, a0, α0, u, v, ω)
    N_v, A_v, M_v = quasisteady1_loads_v(a, b, ρ, a0, α0, u, v, ω)
    N_U = N_u * u_U + N_v * v_U# + L_ω * ω_U
    M_U = M_u * u_U + M_v * v_U# + M_ω * ω_U

    N_ρ, A_ρ, M_ρ = quasisteady1_loads_ρ(a, b, ρ, a0, α0, u, v, ω)

    # return jacobian
    return @SMatrix [
        N_a N_b N_a0 N_α0 N_kh N_kθ N_m N_Sθ N_Iθ N_U N_ρ;
        M_a M_b M_a0 M_α0 M_kh M_kθ M_m M_Sθ M_Iθ M_U M_ρ
        ]
end

function get_coupling_rate_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, p)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_section_rate_jacobian(a, b, ρ)
end

function get_coupling_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, p)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_section_state_jacobian(a, b, ρ, a0, U)
end

function get_coupling_parameter_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, dx, x, p, t)
    # extract state variables
    dh, dθ, dhdot, dθdot = dx
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    u_U, v_U, ω_U = section_velocities_U(θ)
    # local freestream accelerations
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)

    # calculate loads
    N_a, A_a, M_a = quasisteady2_loads_a(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    N_b, A_b, M_b = quasisteady2_loads_b(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    N_a0, A_a0, M_a0 = quasisteady2_loads_a0(a, b, ρ, a0, α0, u, v, ω)
    N_α0, A_α0, M_α0 = quasisteady2_loads_α0(a, b, ρ, a0, α0, u, v, ω)

    N_kh, A_kh, M_kh = 0, 0, 0
    N_kθ, A_kθ, M_kθ = 0, 0, 0
    N_m, A_m, M_m = 0, 0, 0
    N_Sθ, A_Sθ, M_Sθ = 0, 0, 0
    N_Iθ, A_Iθ, M_Iθ = 0, 0, 0

    N_u, A_u, M_u = quasisteady2_loads_u(a, b, ρ, a0, α0, u, v, ω)
    N_v, A_v, M_v = quasisteady2_loads_v(a, b, ρ, a0, α0, u, v, ω)
    N_U = N_u * u_U + N_v * v_U# + L_ω * ω_U
    M_U = M_u * u_U + M_v * v_U# + M_ω * ω_U

    N_ρ, A_ρ, M_ρ = quasisteady2_loads_ρ(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)

    # return jacobian
    return @SMatrix [
        N_a N_b N_a0 N_α0 N_kh N_kθ N_m N_Sθ N_Iθ N_U N_ρ;
        M_a M_b M_a0 M_α0 M_kh M_kθ M_m M_Sθ M_Iθ M_U M_ρ
        ]
end

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- Plotting --- #

@recipe function f(aero::QuasiSteady, stru::TypicalSection, dx, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end

# --- Internal Methods --- #

function section_steady_velocities(U, θ)

    u = U
    v = U*θ

    return SVector(u, v)
end

function section_steady_velocities_U(θ)

    u_U = 1
    v_U = θ

    return SVector(u_U, v_U)
end

function section_steady_velocities_θ(U)

    u_θ = 0
    v_θ = U

    return SVector(u_θ, v_θ)
end


function quasisteady0_section_state_jacobian(a, b, ρ, a0, U)
    L_θ = a0*ρ*U^2*b
    M_θ = (b/2 + a*b)*L_θ
    return @SMatrix [0 L_θ 0 0; 0 M_θ 0 0]
end

function quasisteady1_section_state_jacobian(a, b, ρ, a0, U)
    tmp1 = a0*ρ*U*b
    tmp2 = pi*ρ*b^3
    d1 = b/2 - a*b
    d2 = b/2 + a*b
    L_θ = tmp1*U
    L_hdot = tmp1
    L_θdot = tmp1*d1 + tmp2*U/b
    M_θ = d2*L_θ
    M_hdot = d2*L_hdot
    M_θdot = -tmp2*U + d2*L_θdot
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
end

quasisteady2_section_state_jacobian(a, b, ρ, a0, U) = quasisteady1_section_state_jacobian(a, b, ρ, a0, U)

function quasisteady2_section_rate_jacobian(a, b, ρ)
    N_hddot, A_hddot, M_hddot = quasisteady2_loads_vdot(a, b, ρ)
    N_θddot, A_θddot, M_θddot = quasisteady2_loads_ωdot(a, b, ρ)
    return @SMatrix [0 0 N_hddot N_θddot; 0 0 M_hddot M_θddot]
end
