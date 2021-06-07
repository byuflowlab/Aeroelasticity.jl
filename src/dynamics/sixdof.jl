"""
    RigidBody <: AbstractModel

Six-degree of freedom rigid body dynamic model with state variables ``q =
\\begin{bmatrix} x & y & z & \\phi & \\theta & \\psi & u & v & w & p & q & r
\\end{bmatrix}^T``, inputs ``r = \\begin{bmatrix} F_x & F_y & F_z & M_x & M_y &
M_z \\end{bmatrix}^T`` and parameters ``p = \\begin{bmatrix} m & I_{xx} & I_{yy}
& I_{zz} & I_{xz} & I_{xy} & I_{yz} \\end{bmatrix}^T``.
"""
struct RigidBody <: AbstractModel end

# --- Traits --- #
number_of_states(::Type{RigidBody}) = 12
number_of_inputs(::Type{RigidBody}) = 6
number_of_parameters(::Type{RigidBody}) = 7
inplaceness(::Type{RigidBody}) = OutOfPlace()
mass_matrix_type(::Type{RigidBody}) = Identity()
state_jacobian_type(::Type{RigidBody}) = Varying()
input_jacobian_type(::Type{RigidBody}) = Varying()
input_dependence_type(::Type{RigidBody}) = Linear()

# --- Methods --- #

function get_rates(::RigidBody, q, r, p, t)
    # extract states
    x, y, z, phi, theta, psi, u, v, w, p, q, r = q
    # extract inputs
    Fx, Fy, Fz, Mx, My, Mz = r
    # extract parameters
    m, I_xx, I_yy, I_zz, I_xz, I_xy, I_yz = p
    # calculate state rates
    return rigid_rates()
end

function get_input_jacobian(::RigidBody, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return rigid_input_jacobian(mass, I_xx, I_yy, I_zz, I_xz, I_xy, I_yz)
end

# TODO: Add parameter jacobian

# --- Internal --- #

function rigid_rates(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r, Fx, Fy, Fz, Mx, My, Mz,
    mass, I_xx, I_yy, I_zz, I_xz, I_xy, I_yz)

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
    vdot = F/mass - cross(Ωb, Vb)

    # angular dynamics
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot = I \ (M - cross(Ωb, I*Ωb))

    return SVector(rdot..., ϕdot, θdot, ψdot, vdot..., ωdot...)
end

function rigid_input_jacobian(mass, I_xx, I_yy, I_zz, I_xz, I_xy, I_yz)

    rdot_F = @SMatrix zeros(3, 3)
    ϕdot_F = @SMatrix zeros(1, 3)
    θdot_F = @SMatrix zeros(1, 3)
    ψdot_F = @SMatrix zeros(1, 3)
    vdot_F = (@SMatrix [1 0 0; 0 1 0; 0 0 1])/mass
    ωdot_F = @SMatrix zeros(3, 3)

    rdot_M = @SMatrix zeros(3, 3)
    ϕdot_M = @SMatrix zeros(1, 3)
    θdot_M = @SMatrix zeros(1, 3)
    ψdot_M = @SMatrix zeros(1, 3)
    vdot_M = @SMatrix zeros(3, 3)
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot_M = I \ (@SMatrix [1 0 0; 0 1 0; 0 0 1])

    return [rdot_F rdot_M; ϕdot_F ϕdot_M; θdot_F, θdot_M; ψdot_F, ψdot_M,
        vdot_F vdot_M; ωdot_F ωdot_M]
end
