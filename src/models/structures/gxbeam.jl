"""
    GXBeamAssembly{S}

Model which describes the dynamic behavior of an assembly of beam elements, as modeled by 
the `GXBeam` package.  State variables are as defined by GXBeam.  Inputs are defined and
passed as a [`GXBeamInputs`](@ref) struct. Parameters are defined and passed as a
[`GXBeamParameters`](@ref) struct.  
"""
struct GXBeamAssembly{S}
    system::S
    steady_state::Bool
    two_dimensional::Bool
    structural_damping::Bool
end

"""
    GXBeamAssembly(system::GXBeam.AbstractSystem; steady_state=false, two_dimensional=false, 
        structural_damping=true)

Construct a geometrically exact beam theory structural model
"""
function GXBeamAssembly(system::GXBeam.AbstractSystem;
    steady_state=false, two_dimensional=false, structural_damping=true)

    GXBeamAssembly(system, steady_state, two_dimensional, structural_damping)
end

# residual function
function (gxbeam::GXBeamAssembly)(resid, dx, x, y, p, t)

    # extract constants
    system = gxbeam.system
    indices = system.indices
    force_scaling = system.force_scaling
    steady_state = gxbeam.steady_state
    two_dimensional = gxbeam.two_dimensional
    structural_damping = gxbeam.structural_damping

    # extract parameters
    assembly = p

    # extract inputs
    prescribed, distributed, point_masses, linear_velocity, angular_velocity, gravity = y

    # extract additional inputs for steady simulations
    linear_acceleration, angular_acceleration = y

    # compute residual
    if typeof(system) <: StaticSystem
        static_system_residual!(resid, x, indices, two_dimensional, force_scaling,
            assembly, prescribed, distributed, point_masses, gravity)
    elseif typeof(system) <: DynamicSystem
        if steady_state
            steady_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
                structural_damping, assembly, prescribed, distributed, point_masses,
                gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
        else
            dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
                structural_damping, assembly, prescribed, distributed, point_masses,
                gravity, linear_velocity, angular_velocity)
        end
    else # typeof(system) <: ExpandedSystem
        if steady_state
            expanded_steady_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
                structural_damping, assembly, prescribed, distributed, point_masses,
                gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
        else
            expanded_dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
                structural_damping, assembly, prescribed, distributed, point_masses,
                gravity, linear_velocity, angular_velocity)
        end
    end

    return resid
end

number_of_states(gxbeam::GXBeamAssembly) = length(gxbeam.system.x)

struct GXBeamInputs{TF}
    prescribed_conditions::Dict{Int,PrescribedConditions{TF}}
    distributed_loads::Dict{Int,DistributedLoads{TF}}
    point_masses::Dict{Int,PointMass{TF}}
    linear_velocity::SVector{3,TF}
    angular_velocity::SVector{3,TF}
    linear_acceleration::SVector{3,TF}
    angular_acceleration::SVector{3,TF}
    gravity::SVector{3,TF}
end

"""
    GXBeamInputs(; kwargs...)

Defines inputs for a geometrically exact beam theory structural model

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
        A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and values of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
        with keys corresponding to the elements to which distributed loads are
        applied and values of type [`DistributedLoads`](@ref) which describe
        the distributed loads on those elements.  
 - `point_masses = Dict{Int,PointMass{Float64}}()`: A dictionary with keys
        corresponding to the points to which point masses are attached and values
        of type [`PointMass`](@ref) which contain the properties of the attached
        point masses.
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the body frame.
"""
function GXBeamInputs(;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)))

    return GXBeamInputs(prescribed_conditions, distributed_loads, point_masses, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
        gravity)
end

struct GXBeamParameters{TF,TP,TC,TE}
    assembly::Assembly{TF,TP,TC,TE}
end

"""
    GXBeamParameters(assembly::Assembly)

Defines parameters for a geometrically exact beam theory structural model
"""
GXBeamParameters(assembly)

# --- Internal Methods for Couplings with this Model --- #

# transformation from body to local frame
gxbeam_body_to_element(elem, θ) = elem.Cab'*GXBeam.get_C(θ)

# local element velocities (in local frame)
function gxbeam_element_velocities(V, Ω)
    vi = V
    ωi = Ω
    return vi, ωi
end

# local element accelerations (in local frame)
function gxbeam_element_accelerations(elem, u, θ, Vdot, Ωdot, a0, α0)
    x0 = @SVector zeros(3)
    C = gxbeam_body_to_element(elem, θ)
    ai = C*(a0 + cross(α0, elem.x - x0) + cross(α0, u)) + Vdot
    αi = C*α0 + Ωdot
    return ai, αi
end