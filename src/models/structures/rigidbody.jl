"""
    RigidBody(state_indices=(), rate_indices=(), prescribed_values=())

Construct a six-degree of freedom rigid-body model with state variables ``xr, yr, zr, \\phi,
\\theta, \\psi, ur, vr, wr, pr, qr, rr`` and inputs ``m, I_{xx}, I_{yy}, I_{zz}, I_{xz},
I_{xy}, I_{yz}, F_x, F_y, F_z, M_x, M_y, M_z``.  Uses the state 
variables corresponding to `state_indices` to set the state rates corresponding to
`rate_indices` to the values specified in `prescribed_values`. Otherwise, allow the state 
variables and their respective rates to be defined by their rate equations.
"""
struct RigidBody{IS, IR, V}
    state_indices::IS
    rate_indices::IR
    prescribed_values::V

    function RigidBody(state_indices=(), rate_indices=(), prescribed_values=())
        @assert length(state_indices) == length(rate_indices) == length(prescribed_values)
        IS = typeof(state_indices)
        IR = typeof(rate_indices)
        V = typeof(prescribed_values)
        return new{IS,IR,V}(state_indices, rate_indices, prescribed_values)
    end

end

# residual function
function (rigid_body::RigidBody)(resid, dx, x, y, p, t)
    
    # unpack constants
    @unpack state_indices, rate_indices, prescribed_values = rigid_body

    # extract states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = x
    
    # extract inputs
    m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = y
    
    # rigid body kinematics
    xrdot, yrdot, zrdot, ϕrdot, θrdot, ψrdot = rigidbody_kinematics(xr, yr, zr, ϕr, θr, ψr, 
        ur, vr, wr, pr, qr, rr)
    
    # rigid body dynamics
    urdot, vrdot, wrdot, prdot, qrdot, rrdot = rigidbody_dynamics(ur, vr, wr, pr, qr, rr, m,
        Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)
    
    # define default residuals
    resid[1] = dx[1] - xrdot
    resid[2] = dx[2] - yrdot
    resid[3] = dx[3] - zrdot
    resid[4] = dx[4] - ϕrdot
    resid[5] = dx[5] - θrdot
    resid[6] = dx[6] - ψrdot
    resid[7] = dx[7] - urdot
    resid[8] = dx[8] - vrdot
    resid[9] = dx[9] - wrdot
    resid[10] = dx[10] - prdot
    resid[11] = dx[11] - qrdot
    resid[12] = dx[12] - rrdot

    # replace elements of the residual vector with new constraints
    for (istate, irate, value) in zip(state_indices, rate_indices, prescribed_values)
        residual[istate] = value - rates[irate]
    end
    
    # return result
    return residual
end

# rigid body kinematics
function rigidbody_kinematics(xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr)

    Vb = SVector(ur, vr, wr)
    Ωb = SVector(pr, qr, rr)

    sϕ, cϕ = sincos(ϕr)
    sθ, cθ = sincos(θr)
    sψ, cψ = sincos(ψr)

    # linear kinematics
    Rib = @SMatrix [cθ*cψ    cθ*sψ           -sθ;
         sϕ*sθ*cψ - cϕ*sψ sϕ*sθ*sψ + cϕ*cψ sϕ*cθ;
         cϕ*sθ*cψ + sϕ*sψ cϕ*sθ*sψ - sϕ*cψ cϕ*cθ]
    rdot = Rib' * Vb

    # angular kinematics
    ϕdot = pr + (qr*sϕ + rr*cϕ)*sθ/cθ
    θdot = qr*cϕ - rr*sϕ
    ψdot = (qr*sϕ + rr*cϕ)/cθ

    return SVector(rdot..., ϕdot, θdot, ψdot)
end

# rigid body dynamics
function rigidbody_dynamics(ur, vr, wr, pr, qr, rr, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz,
    Fx, Fy, Fz, Mx, My, Mz)

    F = SVector(Fx, Fy, Fz)
    M = SVector(Mx, My, Mz)

    Vb = SVector(ur, vr, wr)
    Ωb = SVector(pr, qr, rr)

    # linear dynamics
    vdot = F/m - cross(Ωb, Vb)

    # angular dynamics
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot = I \ (M - cross(Ωb, I*Ωb))

    return SVector(vdot..., ωdot...)
end

# --- Internal Methods for Couplings with this Model --- #

function rigidbody_velocities(rigidbody_states)

    v = SVector(rigidbody_states.ur, rigidbody_states.vr, rigidbody_states.wr)
    ω = SVector(rigidbody_states.pr, rigidbody_states.qr, rigidbody_states.rr)
    
    return v, ω
end

function rigidbody_accelerations(rigidbody_rates)

    a = SVector(rigidbody_rates.ur, rigidbody_rates.vr, rigidbody_rates.wr)
    alpha = SVector(rigidbody_rates.pr, rigidbody_rates.qr, rigidbody_rates.rr)

    return a, alpha
end