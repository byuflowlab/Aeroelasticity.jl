"""
    LiftingLineGXBeamAssembly()

Construct a model by coupling a lifting line aerodynamic model and a geometrically exact
beam theory model (as implemented by GXBeam).  Parameters for this model are defined in
[`LiftingLineGXBeamParameters`](@ref).

**NOTE: When using this model, the local frame for each beam element is assumed to
be oriented with the x-axis along the beam's axis, the y-axis forward (into the freestream),
and the z-axis in the airfoil normal direction.**
"""
struct LiftingLineGXBeamAssembly{R, Y, I, S}
    liftingline::LiftingLine{R, Y, I}
    gxbeam::GXBeamAssembly{S}
end

# set as default coupling
default_coupling(liftingline::LiftingLine, gxbeam::GXBeamAssembly) = LiftingLineGXBeamAssembly(liftingline, gxbeam)

# default coupling function
function (liftingline_gxbeam::LiftingLineGXBeamAssembly)(dx, x, p, t)

    # coupling function for each section
    finput = liftingline_gxbeam.liftingline.finput

    # state variable indices for each section
    indices = liftingline_gxbeam.liftingline.indices

    # extract rates
    liftingline_rates, gxbeam_rates = dx

    # extract states
    liftingline_states, gxbeam_states = x

    # extract parameters
    liftingline_parameters, gxbeam_parameters, coupling_parameters = p

    # extract GXBeam variables
    system = liftingline_gxbeam.gxbeam.system
    assembly = gxbeam_parameters
    prescribed_conditions = coupling_parameters.prescribed_conditions
    distributed_loads = coupling_parameters.distributed_loads

    # number of sections
    n = length(indices)

    # freestream properties
    Vinf = coupling_parameters.Vinf
    rho = coupling_parameters.rho
    beta = coupling_parameters.beta

    # initialize aerodynamic model inputs
    liftingline_inputs = ()

    # initialize combined distributed loads
    combined_distributed_loads = ()

    # loop through each section
    for i = 1:n

        # extract element properties
        element = assembly.elements[i]

        # extract element state variables
        element_states = GXBeam.extract_element_state(gxbeam_rates, gxbeam_states,
            system, assembly, i; prescribed_conditions = prescribed_conditions)

        # construct transformation matrices for this element
        body_to_element = gxbeam_body_to_element(element, element_states.theta)
        element_to_aero = @SMatrix [0 -1 0; 1 0 0; 0 0 1]
        body_to_aero = element_to_aero*body_to_element

        # local freestream velocity and acceleration
        vi = body_to_aero * (Vinf - element_states.V)
        ωi = -body_to_aero*element_states.Omega
        ai = -body_to_aero*element_states.Vdot
        αi = -body_to_aero*element_states.Omegadot

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

        # # transform distributed loads into the body frame
        fi = body_to_aero'*fi
        mi = body_to_aero'*mi

        # integrate distributed loads
        section_loads = DistributedLoads(assembly, i;
            fx = (s) -> fi[1], fy = (s) -> fi[2], fz = (s) -> fi[3],
            mx = (s) -> mi[1], my = (s) -> mi[2], mz = (s) -> mi[3])

        # combine with existing distributed loads
        if haskey(distributed_loads, i)
            section_loads = GXBeam.combine_loads(distributed_loads[i], section_loads)
        end

        # save distributed loads
        combined_distributed_loads = (combined_distributed_loads..., section_loads)
    end

    # construct gxbeam inputs
    gxbeam_inputs = GXBeamInputs(;
        prescribed_conditions = coupling_parameters.prescribed_conditions,
        distributed_loads = Dict(i=>combined_distributed_loads[i] for i = 1:n),
        point_masses = coupling_parameters.point_masses,
        linear_velocity = coupling_parameters.linear_velocity,
        angular_velocity = coupling_parameters.angular_velocity,
        linear_acceleration = coupling_parameters.linear_acceleration,
        angular_acceleration = coupling_parameters.angular_acceleration,
        gravity = coupling_parameters.gravity)

    return liftingline_inputs, gxbeam_inputs
end

"""
    LiftingLineGXBeamParameters(Vinf, rho, beta; kwargs...)

Defines parameters for a lifting line model coupled with a geometrically exact beam theory
structural model when subjected to a the freestream velocity vector `Vinf`.

# Arguments
 - `Vinf`: Freestream velocity vector
 - `rho`: Air density
 - `beta`: Prandtl-Glauert compressibility factor ``\\beta = \\sqrt{1 - M^2}``

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
        A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and values of type
        `PrescribedConditions` which describe the prescribed conditions
        at those points.
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
        with keys corresponding to the elements to which distributed loads are
        applied and values of type `DistributedLoads` which describe
        the distributed loads on those elements.
 - `point_masses = Dict{Int,PointMass{Float64}}()`: A dictionary with keys
        corresponding to the points to which point masses are attached and values
        of type `PointMass` which contain the properties of the attached
        point masses.
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the body frame.
"""
function LiftingLineGXBeamParameters(Vinf, rho, beta;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)))

    return (; Vinf, rho, beta, prescribed_conditions, distributed_loads, point_masses,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        gravity)
end