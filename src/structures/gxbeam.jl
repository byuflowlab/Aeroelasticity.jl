"""
    GEBT <: AbstractModel

Geometrically exact beam theory model, as implemented by the GXBeam package.
State variables are as defined by GXBeam. Inputs correspond to the (non-follower)
prescribed forces/displacements on each node with prescribed conditions followed
by the (non-follower) distributed loads on each beam element with distributed loads.
At this point in time, this model doesn't accept any parameters.

When coupled with an aerodynamic model, the local beam y and z-axes should be
aligned with the negative chordwise and positive normal directions, respectively.
"""
struct GEBT{TF, TV, TM} <: AbstractModel
    system::GXBeam.System{TF, TV, TM}
    assembly::GXBeam.Assembly{TF}
    prescribed::Dict{Int,GXBeam.PrescribedConditions{TF}}
    distributed::Dict{Int,GXBeam.DistributedLoads{TF}}
end

# --- Constructors --- #

"""
    GEBT(system, assembly, prescribed_conditions, distributed_loads)

Construct a geometrically exact beam theory structural model.
"""
function GEBT(system::GXBeam.System{TF, TV, TM}, assembly, prescribed,
    distributed) where {TF, TV, TM}

    return GEBT{TF, TV, TM}(system, assembly, prescribed, distributed)
end

# --- Traits --- #

number_of_states(model::GEBT) = length(model.system.x)

function number_of_inputs(model::GEBT)
    return 6*length(model.distributed)
end

number_of_parameters(model::GEBT) = 0

inplaceness(::Type{<:GEBT}) = InPlace()

mass_matrix_type(::Type{<:GEBT}) = Linear()

state_jacobian_type(::Type{<:GEBT}) = Nonlinear()

input_jacobian_type(::Type{<:GEBT}) = Linear()

# --- Methods --- #

function get_rates!(dq, model::GEBT, q, r, p, t)

    # extract system, assembly, prescribed conditions and distributed loads
    system = model.system
    assembly = model.assembly
    prescribed = model.prescribed
    distributed = model.distributed

    # extract system constants and pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # set origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector zeros(3)
    ω0 = @SVector zeros(3)

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(r)
        # use pre-allocated storage
        update_distributed!(distributed, assembly.elements, r)
    else
        # allocate new storage
        distributed = update_distributed(distributed, assembly.elements, r)
    end

    # return mass matrix multiplied state rates
    return gxbeam_rates!(dq, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_mass_matrix!(M, model::GEBT, q, r, p, t)

    # extract system and assembly
    system = model.system
    assembly = model.assembly

    # extract system constants and  pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # return mass matrix
    return gxbeam_mass_matrix!(M, q, assembly, force_scaling, mass_scaling,
        irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
end

# --- Performance Overloads --- #

function get_state_jacobian!(J, model::GEBT, q, r, p, t)

    # extract system, assembly, prescribed conditions and distributed loads
    system = model.system
    assembly = model.assembly
    prescribed = model.prescribed
    distributed = model.distributed

    # extract system constants and  pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # set origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector zeros(3)
    ω0 = @SVector zeros(3)

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(r)
        # use pre-allocated storage
        update_distributed!(distributed, assembly.elements, r)
    else
        # allocate new storage
        distributed = update_distributed(distributed, assembly.elements, r)
    end

    # return jacobian of right hand side with respect to state variables
    return gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_input_jacobian(Jr, model::GEBT, q, r, p, t)

    # extract system, assembly, and distributed loads
    system = model.system
    assembly = model.assembly
    distributed = model.distributed

    # extract start and stop of each beam element
    elements = assembly.elements
    start = assembly.start
    stop = assembly.stop

    # extract system constants and pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt

    # get elements to which distributed loads are applied
    element_indices = keys(distributed)

    # return input jacobian linear map
    return gxbeam_input_jacobian(q, r, elements, start, stop, force_scaling,
        mass_scaling, irow_pt, element_indices)
end

# --- Unit Testing Methods --- #

function get_lhs(model::GEBT, dq, q, r, p, t)

    # extract system, assembly, prescribed conditions and distributed loads
    system = model.system
    assembly = model.assembly
    prescribed = model.prescribed
    distributed = model.distributed

    # extract system constants and  pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # set origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector zeros(3)
    ω0 = @SVector zeros(3)

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(r)
        # use pre-allocated storage
        update_distributed!(distributed, assembly.elements, r)
    else
        # allocate new storage
        distributed = update_distributed(distributed, assembly.elements, r)
    end

    return gxbeam_lhs(q, dq, assembly, prescribed, distributed, force_scaling,
        mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, x0, v0, ω0)
end

# --- Internal --- #

function gxbeam_lhs(q, dq, assembly, prescribed, distributed,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2,
    icol_pt, icol_beam, x0, v0, ω0)

    steady_residual = similar(dq)
    dynamic_residual = similar(dq)

    GXBeam.steady_state_system_residual!(steady_residual, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    GXBeam.dynamic_system_residual!(dynamic_residual, q, dq, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return steady_residual - dynamic_residual
end

function gxbeam_rhs!(out, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    # mass matrix multiplied rates are equal to steady state GXBeam residuals
    GXBeam.steady_state_system_residual!(out, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return out
end

function gxbeam_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
    irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    M .= 0

    # mass matrix is GXBeam mass matrix, moved to LHS
    GXBeam.system_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    M .*= -1

    return M
end

function gxbeam_rates!(Mdu, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    gxbeam_rhs!(Mdu, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return Mdu
end

function gxbeam_state_jacobian!(K, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    K .= 0

    # jacobian of mass matrix multiplied rates is equal to steady state jacobian
    GXBeam.steady_state_system_jacobian!(K, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return K
end

function gxbeam_input_jacobian(q, r, elements, start, stop, force_scaling, mass_scaling,
    irow_pt, element_indices)

    f! = (y, x) -> gxbeam_jac_vec!(y, x, r, elements, start, stop, force_scaling,
        mass_scaling, irow_pt, element_indices)
    M = length(q)
    N = length(r)

    return LinearMap(f!, M, N; ismutating = true)
end

# input jacobian vector product
function gxbeam_jac_vec!(y, x, r, elements, start, stop, force_scaling, mass_scaling,
    irow_pt, element_indices)

    nelem = length(start)

    # initialize output array
    y .= 0

    # current input array index
    ir = 0

    # distributed loads
    for ielem = 1:nelem
        if ielem in element_indices
            # get point indices
            irow_p1 = irow_pt[start[ielem]]
            irow_p2 = irow_pt[stop[ielem]]
            # initialize element load jacobians
            f_r = Diagonal((@SVector ones(3)))
            m_r = Diagonal((@SVector ones(3)))
            # jacobians after integration
            # (assume distributed loads are constant on each element)
            f1_r = f2_r = elements[ielem].L*f_r/2
            m1_r = m2_r = elements[ielem].L*m_r/2
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
function update_distributed(distributed, elements, r::AbstractVector{T}) where T
    new_distributed = Dict{Int, DistributedLoads{T}}()
    for (key, value) in distributed
        new_distributed[key] = DistributedLoads{T}(value)
    end
    update_distributed!(new_distributed, elements, r)
    return new_distributed
end

# update distributed loads using inputs
function update_distributed!(distributed, elements, r)
    distributed_points = keys(distributed)
    largest_point = maximum(distributed_points)
    ir = 0
    for ielem in 1:largest_point
        if ielem in distributed_points
            # extract loads on this element
            f = SVector(r[ir+1], r[ir+2], r[ir+3])
            m = SVector(r[ir+4], r[ir+5], r[ir+6])
            # get integrated distributed loads
            f1 = f2 = elements[ielem].L*f/2
            m1 = m2 = elements[ielem].L*m/2
            # keep the same follower loads
            f1_follower = distributed[ielem].f1_follower
            f2_follower = distributed[ielem].f2_follower
            m1_follower = distributed[ielem].m1_follower
            m2_follower = distributed[ielem].m2_follower
            # update distributed loads
            distributed[ielem] = DistributedLoads(f1, f2, m1, m2, f1_follower,
                f2_follower, m1_follower, m2_follower)
            ir += 6
        end
    end

    return distributed
end
