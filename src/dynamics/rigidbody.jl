"""
    RigidBody <: AbstractModel

Six-degree of freedom rigid-body model with state variables ``x, y, z, \\phi,
\\theta, \\psi, u, v, w, p, q, r``, inputs ``m, I_{xx}, I_{yy}, I_{zz}, I_{xz},
I_{xy}, I_{yz}, F_x, F_y, F_z, M_x, M_y, M_z``, and zero parameters.
"""
struct RigidBody{N,TI,TF} <: AbstractModel
    state_indices::NTuple{N,TI}
    rate_indices::NTuple{N,TI}
    rate_values::NTuple{N,TF}
end

# --- Constructors --- #

"""
    RigidBody()

Initialize an object of type [`RigidBody`](@ref)
"""
RigidBody() = RigidBody{0,Int,Float64}((), (), ())

"""
    RigidBody(state_indices, rate_indices [, rate_values])

Initialize an object of type [`RigidBody`](@ref).  Use the state variables
corresponding to `state_indices` to set the state rates corresponding to
`rate_indices` to the values specified in `rate_values` (which default to zeros
if not specified). Otherwise, allow the state variables and their respective
rates to be defined by their rate equations.
"""
function RigidBody(state_indices::NTuple{N,TI}, rate_indices::NTuple{N,TI}) where {N,TI}

    rate_values = ntuple(i->0.0, N)

    return RigidBody(state_indices, rate_indices, rate_values)
end

# --- Traits --- #
number_of_states(::Type{<:RigidBody}) = 12
number_of_inputs(::Type{<:RigidBody}) = 13
number_of_parameters(::Type{<:RigidBody}) = 0
inplaceness(::Type{<:RigidBody}) = OutOfPlace()
mass_matrix_type(::Type{RigidBody{0,TI,TF}}) where {TI,TF} = Identity()
mass_matrix_type(::Type{RigidBody{N,TI,TF}}) where {N,TI,TF} = Linear()
state_jacobian_type(::Type{<:RigidBody}) = Nonlinear()
input_jacobian_type(::Type{<:RigidBody}) = Nonlinear()

# --- Methods --- #

function get_rates(model::RigidBody, states, inputs, parameters, t)
    # extract states
    x, y, z, ϕ, θ, ψ, u, v, w, p, q, r = states
    # extract inputs
    m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = inputs
    # calculate state rates
    dx, dy, dz, dϕ, dθ, dψ = rigid_body_kinematics(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r)
    du, dv, dw, dp, dq, dr = rigid_body_dynamics(u, v, w, p, q, r, m,
        Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)
    # impose constraints
    rates = SVector(dx, dy, dz, dϕ, dθ, dψ, du, dv, dw, dp, dq, dr)
    for (istate, irate, value) in zip(model.state_indices, model.rate_indices, model.rate_values)
        rates = setindex(rates, value, istate)
    end
    # return result
    return rates
end

function get_mass_matrix(model::RigidBody, states, inputs, parameters, t)
    # start with state rate equation mass matrix
    M = SMatrix{12, 12, TF}(I)
    # impose constraints
    for (istate, irate, value) in zip(model.state_indices, model.rate_indices, model.rate_values)
        M = setindex(M, value, istate, irate)
    end
    # return result
    return M
end

# --- Performance Overloads --- #

#TODO: Add state jacobian

#TODO: Add input jacobian

# --- Unit Testing Methods --- #

function get_lhs(model::RigidBody, dstates, states, inputs, parameters, t)
    # impose constraints
    dstates = SVector{12}(dstates)
    for (istate, irate, value) in zip(model.state_indices, model.rate_indices, model.rate_values)
        dstates = setindex(dstates, dstates[irate], istate)
    end
    return dstates
end

# --- Internal Methods --- #

function rigid_body_kinematics(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r)

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

    return SVector(rdot..., ϕdot, θdot, ψdot)
end

function rigid_body_dynamics(u, v, w, p, q, r, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz,
    Fx, Fy, Fz, Mx, My, Mz)

    F = SVector(Fx, Fy, Fz)
    M = SVector(Mx, My, Mz)

    Vb = SVector(u, v, w)
    Ωb = SVector(p, q, r)

    # linear dynamics
    vdot = F/m - cross(Ωb, Vb)

    # angular dynamics
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot = I \ (M - cross(Ωb, I*Ωb))

    return SVector(vdot..., ωdot...)
end
