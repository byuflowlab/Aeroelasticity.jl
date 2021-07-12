"""
    RigidBody <: AbstractModel

Six-degree of freedom rigid body dynamic model with state variables ``q =
\\begin{bmatrix} x & y & z & \\phi & \\theta & \\psi & u & v & w & p & q & r
\\end{bmatrix}^T``, inputs ``r = \\begin{bmatrix} m & I_{xx} & I_{yy}
& I_{zz} & I_{xz} & I_{xy} & I_{yz} & F_x & F_y & F_z & M_x & M_y &
M_z \\end{bmatrix}^T`` and zero parameters.
"""
struct RigidBody <: AbstractModel end

"""
    Gravity <: NoStateModel

Gravitational model with parameter ``g``
"""
struct Gravity <: NoStateModel end

# --- Traits --- #
number_of_states(::Type{RigidBody}) = 12
number_of_inputs(::Type{RigidBody}) = 13
number_of_parameters(::Type{RigidBody}) = 0
inplaceness(::Type{RigidBody}) = OutOfPlace()
mass_matrix_type(::Type{RigidBody}) = Identity()
state_jacobian_type(::Type{RigidBody}) = Nonlinear()
input_jacobian_type(::Type{RigidBody}) = Nonlinear()

number_of_parameters(::Type{Gravity}) = 1

# --- Methods --- #

function get_rates(::RigidBody, states, inputs, parameters, t)
    # extract states
    x, y, z, ϕ, θ, ψ, u, v, w, p, q, r = states
    # extract inputs
    m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = inputs
    # calculate state rates
    return rigid_rates(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r, m, Ixx, Iyy, Izz,
        Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)
end

# --- Performance Overloads --- #

#TODO: Add state jacobian

#TODO: Add input jacobian

# --- Unit Testing Methods --- #

get_lhs(::RigidBody, dstates, states, inputs, parameters, t) = dstates

# --- Internal Methods --- #

function rigid_rates(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r, m, Ixx, Iyy, Izz, Ixz,
    Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)

    F = SVector(Fx, Fy, Fz)
    M = SVector(Mx, My, Mz)

    Vb = SVector(u, v, w)
    Ωb = SVector(p, q, r)

    sϕ, cϕ = sincos(ϕ)
    sθ, cθ = sincos(θ)
    sψ, cψ = sincos(ψ)

    # linear kinematics
    Rib = @SMatrix [cθ*cψ    cθ*sψ           -sθ;
         sϕ*sθ*cψ - cϕ*sψ sϕ*sθ*sψ + cϕ*cψ sϕ*cθ;
         cϕ*sθ*cψ + sϕ*sψ cϕ*sθ*sψ - sϕ*cψ cϕ*cθ]
    rdot = Rib' * Vb

    # angular kinematics
    ϕdot = p + (q*sϕ + r*cϕ)*sθ/cθ
    θdot = q*cϕ - r*sϕ
    ψdot = (q*sϕ + r*cϕ)/cθ

    # linear dynamics
    vdot = F/m - cross(Ωb, Vb)

    # angular dynamics
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot = I \ (M - cross(Ωb, I*Ωb))

    return SVector(rdot..., ϕdot, θdot, ψdot, vdot..., ωdot...)
end
