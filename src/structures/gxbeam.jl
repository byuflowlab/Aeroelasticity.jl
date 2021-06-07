"""
    GeometricallyExactBeamTheory <: AbstractModel

Geometrically exact beam theory model, as implemented by the GXBeam package.
State variables are as defined by GXBeam. Inputs correspond to the (non-follower)
distributed aerodynamic loads on each beam element with distributed loads.
Currently, this model doesn't accept any parameters.
"""
struct GeometricallyExactBeamTheory{TF, TV, TM} <: AbstractModel
    system::GXBeam.System{TF, TV, TM}
    assembly::GXBeam.Assembly{TF}
    prescribed_conditions::Dict{Int,GXBeam.PrescribedConditions{TF}}
    distributed_loads::Dict{Int,GXBeam.DistributedLoads{TF}}
end

# --- Constructors --- #

"""
    GeometricallyExactBeamTheory(system, assembly, prescribed_conditions, distributed_loads)

Construct a GeometricallyExactBeamTheory structural model.
"""
GeometricallyExactBeamTheory(system, assembly, prescribed_conditions, distributed_loads)

# --- Traits --- #

number_of_states(model::GeometricallyExactBeamTheory) = length(model.system.x)
function number_of_inputs(model::GeometricallyExactBeamTheory)
    return 6*length(model.distributed)
end
number_of_parameters(model::GeometricallyExactBeamTheory) = 0
inplaceness(::Type{GeometricallyExactBeamTheory}) = InPlace()
mass_matrix_type(::Type{GeometricallyExactBeamTheory}) = Varying()
state_jacobian_type(::Type{GeometricallyExactBeamTheory}) = Varying()
input_jacobian_type(::Type{GeometricallyExactBeamTheory}) = Varying()
input_dependence_type(::Type{GeometricallyExactBeamTheory}) = Linear()

# --- Methods --- #

function get_mass_matrix!(M, model::GeometricallyExactBeamTheory, q, r, p, t)

    # extract system pointers
    system = model.system
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # extract system constants
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # extract template assembly
    assembly = model.assembly

    # return mass matrix
    return gebt_mass_matrix!(M, q, assembly, force_scaling, mass_scaling,
        irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
end

function get_rates!(dq, model::GeometricallyExactBeamTheory, q, r, p, t)

    # extract system pointers
    system = model.system
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # extract system constants
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # extract assembly, prescribed conditions and distributed loads
    assembly = model.assembly
    prescribed = model.prescribed
    distributed = model.distributed

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(r)
        # use pre-allocated storage
        update_distributed!(distributed, r)
    else
        # allocate new storage
        distributed = update_distributed(distributed, r)
    end

    # set origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector zeros(3)
    ω0 = @SVector zeros(3)

    # return mass matrix multiplied state rates
    return gxbeam_rates!(Mdu, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_state_jacobian!(J, model::GeometricallyExactBeamTheory, q, r, p, t)

    # extract system pointers
    system = model.system
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # extract system constants
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # extract assembly, prescribed conditions and distributed loads
    assembly = model.assembly
    prescribed = model.prescribed
    distributed = model.distributed

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(r)
        # use pre-allocated storage
        update_distributed!(distributed, r)
    else
        # allocate new storage
        distributed = update_distributed(distributed, r)
    end

    # set origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector zeros(3)
    ω0 = @SVector zeros(3)

    # return jacobian of right hand side with respect to state variables
    return gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_input_jacobian(Jr, model::GeometricallyExactBeamTheory, q, r, p, t)

    # extract system pointers
    system = model.system
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # extract system constants
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # return input jacobian linear map
    return gxbeam_input_jacobian(q, r, start, stop, distributed_loads, irow_pt,
        force_scaling, mass_scaling)
end

# --- Internal --- #

# jacobian of LHS wrt state rates
function gxbeam_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
    irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    system_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    M .*= -1

    return M
end

# jacobian of RHS
function gxbeam_rates!(Mdu, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    GXBeam.steady_state_system_residual!(Mdu, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return Mdu
end

# jacobian of RHS wrt states
function gxbeam_state_jacobian!(K, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    GXBeam.steady_state_system_jacobian!(K, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return K
end

# jacobian of RHS wrt inputs
function gxbeam_input_jacobian(q, r, start, stop, distributed_loads, irow_pt,
    force_scaling, mass_scaling)

    f! = (y, x) -> gxbeam_jac_vec!(y, x, r, start, stop, distributed_loads, irow_pt,
        force_scaling, mass_scaling)
    M = length(q)
    N = length(r)

    return LinearMap(f!, M, N; ismutating = true)
end

# input jacobian vector product
function gxbeam_jac_vec!(y, x, r, start, stop, distributed_loads, irow_pt,
    force_scaling, mass_scaling)

    # initialize output array
    y .= 0

    # current input array index
    ir = 0

    # get all elements with distributed loads
    distributed_load_elements = keys(distributed_loads)

    # distributed loads
    for ielem = 1:nelem
        if ielem in distributed_load_elements
            # get point indices
            irow_p1 = irow_pt[start[ielem]]
            irow_p2 = irow_pt[stop[ielem]]
            # initialize element load jacobians
            f_r = Diagonal((@SVector ones(3)))
            m_r = Diagonal((@SVector ones(3)))
            # jacobians after integration
            f1_r = f2_r = ΔL*f_r/2
            m1_r = m2_r = ΔL*m_r/2
            # calculate jacobian vector product
            xf = SVector(x[ir+1], x[ir+2], x[ir+3])
            xm = SVector(x[ir+4], x[ir+5], x[ir+6])
            y[irow_p1:irow_p1+2] .+= -f1_r * xf ./ force_scaling
            y[irow_p1+3:irow_p1+5] .+= -m1_r * xm ./ force_scaling
            y[irow_p2:irow_p2+2] .+= -f2_r * xf ./ force_scaling
            y[irow_p2+3:irow_p2+5] .+= -m2_r * xm ./ force_scaling
            # move to next set of distributed loads
            ir += 6
        end
    end

    return y
end

# update distributed loads using inputs
function update_distributed(distributed, r::AbstractVector{T}) where T
    new_distributed = Dict{Int, DistributedLoads{T}}()
    for (key, value) in distributed
        new_distributed[key] = value
    end
    update_distributed!(distributed, r)
    return new_distributed
end

# update distributed loads using inputs
function update_distributed!(distributed, r)
    distributed_points = keys(distributed)
    largest_point = maximum(distributed_points)
    ir = 0
    for ipoint in 1:largest_point
        if ipoint in distributed_points
            # extract loads on this element
            f = SVector(r[ir+1], r[ir+2], r[ir+3])
            m = SVector(r[ir+4], r[ir+5], r[ir+6])
            # get integrated distributed loads
            f1 = f2 = ΔL*f/2
            m1 = m2 = ΔL*m/2
            # keep the same follower loads
            f1_follower = distributed.f1_follower
            f2_follower = distributed.f2_follower
            m1_follower = distributed.m1_follower
            m2_follower = distributed.m2_follower
            # update distributed loads
            distributed[ipoint] = DistributedLoads(f1, f2, m1, m2, f1_follower,
                f2_follower, m1_follower, m2_follower)
            ir += 6
        end
    end

    return distributed
end
