"""
    LiftingLineRigidBody()

Construct a model by coupling a lifting line aerodynamic model and a rigid body dynamics 
model. Parameters for this model are defined in [`LiftingLineRigidBodyParameters`](@ref).

**NOTE: When using this model, the local frame for each lifting line element is assumed to
be oriented with the x-axis in the chordwise direction, the y-axis in the spanwise 
direction (e.g. out the right wing), and the z-axis in the airfoil normal direction.**

**NOTE: This model does not define a default parameter function for the coupled model**
"""
struct LiftingLineRigidBody{R, Y, I, IS, IR, V}
    liftingline::LiftingLine{R, Y, I}
    rigidbody::RigidBody{IS, IR, V}
end

# set as default coupling
default_coupling(liftingline::LiftingLine, rigidbody::RigidBody) = LiftingLineRigidBody(liftingline, rigidbody)

# default coupling function
function (liftingline_rigidbody::LiftingLineRigidBody)(dx, x, p, t)

    # coupling function for each section
    finput = liftingline_rigidbody.liftingline.finput

    # state variable indices for each section
    indices = liftingline_rigidbody.liftingline.indices

    # extract rates
    _, _, _, _, _, _, dur, dvr, dwr, dpr, dqr, drr = dx[2]

    # extract states
    _, _, _, _, _, _, ur, vr, wr, pr, qr, rr = x[2]

    # extract parameters
    liftingline_parameters, rigidbody_parameters, coupling_parameters = p

    # number of sections
    n = length(indices)

    # velocity and acceleration vectors
    v = SVector(ur, vr, wr)
    ω = SVector(pr, qr, rr)
    a = SVector(dur, dvr, dwr)
    α = SVector(dpr, dqr, drr)

    # initialize loads with body loads
    F = coupling_parameters.Fb
    M = coupling_parameters.Mb

    # extract freestream density and speed of sound
    rho = coupling_parameters.rho
    beta = coupling_parameters.beta

    # initialize aerodynamic model inputs
    liftingline_inputs = ()

    # loop through each section
    for i = 1:n

        # extract local section properties
        L = coupling_parameters.L[i]
        pe = coupling_parameters.pe[i]
        e1 = coupling_parameters.e1[i]
        e2 = coupling_parameters.e2[i]
        e3 = coupling_parameters.e3[i] 
        
        # calculate local section velocities and accelerations
        vi, ωi, ai, αi = liftingline_section_velocities(pe, e1, e2, e3, v, ω, a, α)

        # define section rates, states, and parameters
        dxi = liftingline_rates[indices[i]], SVector(ai..., αi...)
        xi = liftingline_states[indices[i]], SVector(vi..., ωi...)
        pi = liftingline_parameters[i], (rho, beta)
        ti = t

        # calculate section inputs
        yi_aero, yi_load = finput[i](dxi, xi, pi, ti)

        # save aerodynamic inputs
        liftingline_inputs = (liftingline_inputs..., yi_aero)

        # extract load inputs
        fi, mi = yi_load

        # integrate distributed loads
        Fi = L*fi
        Mi = L*mi

        # add section loads to total loads
        F += Fi
        M += cross(pe, Fi) + Mi
    end

    # save rigid body inputs
    rigidbody_inputs = (coupling_parameters.m, coupling_parameters.I, F, M)

    return liftingline_inputs, rigidbody_inputs
end

"""
    LiftingLineRigidBodyParameters{TF}

Struct containing properties for the [`LiftingLineRigidBody`](@ref) coupling.
"""
struct LiftingLineRigidBodyParameters{TF}
    m::TF
    I::SMatrix{36,6,6,TF}
    L::SVector{3,TF}
    p_e::Vector{SVector{3,TF}}
    e_1::Vector{SVector{3,TF}}
    e_2::Vector{SVector{3,TF}}
    e_3::Vector{SVector{3,TF}}
    rho::TF
    c::TF
    Fb::SVector{3,TF}
    Mb::SVector{3,TF}
end

"""
    LiftingLineRigidBodyParameters(m, I, L, p_e [, e_1, e_2, e_3]; kwargs...)

Constructs an object of type [`LiftingLineRigidBodyParameters`](@ref)

# Arguments:
 - `m`: Body mass
 - `I`: Body inertia matrix
 - `L`: Length of each lifting line element
 - `p_e`: Position of each lifting line element
 - `e_1 = fill([1,0,0], length(L))` Unit vector for the x-axis of each lifting line element
 - `e_2 = fill([0,1,0], length(L))` Unit vector for the y-axis of each lifting line element
 - `e_3 = fill([0,0,1], length(L))` Unit vector for the z-axis of each lifting line element

# Keyword Arguments
 - `rho=1.225`: freestream air density
 - `c=343`: air speed of sound
 - `Fb=[0,0,0]`: applied rigid body force (excluding lifting line forces)
 - `Mb=[0,0,0]`: applied rigid body moment (excluding lifting line moments)
"""
function LiftingLineRigidBodyParameters(m, I, L, p_e, e_1=fill(SVector(1,0,0), length(L)), 
    e_2=fill(SVector(0,1,0), length(L)), e_3=fill(SVector(0,0,1), length(L)); 
    rho = 1.225, c=343, Fb=[0,0,0], Mb=[0,0,0])

    LiftingLineRigidBodyParameters(m, I, L, p_e, e_1, e_2, e_3, rho, c, Fb, Mb)
end