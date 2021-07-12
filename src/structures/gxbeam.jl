"""
    GEBT <: AbstractModel

Geometrically exact beam theory model, as implemented by the GXBeam package.
State variables are as defined by GXBeam.  Inputs correspond to the external
forces ``F_{x,i}``, ``F_{y,i}``, ``F_{z,i}``, ``M_{x,i}``, ``M_{y,i}``,
``M_{z,i}`` or displacements ``u_{x,i}``, ``u_{y,i}``, ``u_{z,i}``,
``\theta_{x,i}``, ``\theta_{y,i}``, ``\theta_{z,i}`` applied to each node,
followed by the distributed loads ``f_{x,i}``, ``f_{y,i}``, ``f_{z,i}``,
``m_{x,i}``, ``m_{y,i}``, ``m_{z,i}`` applied to each beam element, followed by
the linear and angular velocity ``V_x``, ``V_y``, ``V_z``, ``\\Omega_x``,
``\\Omega_y``, ``\\Omega_z`` of the system.  Parameters correspond to the
location ``p_{x}``, ``p_{y}``, ``p_{z}`` of each node followed
by each beam element's properties. Each beam element's properties are defined
by a triad which defines the orientation of the beam element ``e_{1,x}``,
``e_{1,y}``, ``e_{1,z}``, ``e_{2,x}``, ``e_{2,y}``, ``e_{2,z}``, ``e_{3,x}``,
``e_{3,y}``, ``e_{3,z}`` followed by the 21 independent entries of the
compliance matrix ``C_{11}``, ``C_{12}``, ``C_{13}``, ``C_{14}``, ``C_{15}``,
``C_{16}``, ``C_{22}``, ``C_{23}``, ``C_{24}``, ``C_{25}``, ``C_{26}``,
``C_{33}``, ``C_{34}``, ``C_{35}``, ``C_{36}``, ``C_{44}``, ``C_{45}``,
``C_{46}``, ``C_{55}``, ``C_{56}``, ``C_{66}``, followed by the beam element's
inertial properties ``\mu``, ``x_{m,2}``, ``x_{m,3}``, ``i_{22}``, ``i_{33}``,
``i_{23}``.

When coupled with an aerodynamic model, the local beam y and z-axes should be
aligned with the negative chordwise and positive normal directions, respectively.
"""
struct GEBT{TF} <: AbstractModel
    force_scaling::TF
    mass_scaling::TF
    irow_pt::Vector{Int}
    irow_beam::Vector{Int}
    irow_beam1::Vector{Int}
    irow_beam2::Vector{Int}
    icol_pt::Vector{Int}
    icol_beam::Vector{Int}
    element_start::Vector{Int}
    element_stop::Vector{Int}
    isforce::Vector{SVector{6, TF}}
end

# --- Constructors --- #

"""
    GEBT(GXBeam.assembly, prescribed = Dict(Int, <:GXBeam.PrescribedConditions))

Construct a geometrically exact beam theory structural model with connectivity
as specified in ``assembly`` and prescribed displacements (rather than forces)
as specified in ``prescribed``.
"""
function GEBT(assembly, prescribed)

    # get the number of points, elements, and element connections for each point
    npoint = length(assembly.points)
    nelem = length(assembly.elements)
    nconn = GXBeam.point_connections(assembly)

    # set the force scaling based on the average compliance matrix value
    compliance_entries = vcat([vcat(elem.C11..., elem.C12..., elem.C22...)
        for elem in assembly.elements]...)
    compliance_nonzero_indices = findall(xi -> abs(xi) > eps(TF), compliance_entries)
    if isempty(compliance_nonzero_indices)
        force_scaling = 1.0
    else
        nonzero_compliance_entries = compliance_entries[compliance_nonzero_indices]
        force_scaling = nextpow(2.0, length(nonzero_compliance_entries)/
            sum(nonzero_compliance_entries)/100)
    end

    # set the mass scaling based on the average inverse mass matrix value
    minv_entries = vcat([vcat(elem.minv11..., elem.minv12..., elem.minv22...)
        for elem in assembly.elements]...)
    minv_nonzero_indices = findall(xi -> abs(xi) > eps(TF), minv_entries)
    if isempty(minv_nonzero_indices)
        mass_scaling = 1.0
    else
        nonzero_minv_entries = minv_entries[minv_nonzero_indices]
        mass_scaling = nextpow(2.0, length(nonzero_minv_entries)/
            sum(nonzero_minv_entries))
    end

    # construct pointers into the system matrix
    n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam =
        GXBeam.system_indices(assembly, 1:npoint, nconn, false)

    # flags indicating whether point condition is a force or displacement
    isload = fill(@SVector ones(6, Bool), npoint)
    for key in keys(prescribed)
        isload[key] = prescribed[key].force
    end

    return GEBT(force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, assembly.start, assembly.stop, isload)
end

# --- Traits --- #

number_of_states(model::GEBT) = length(model.system.x)

function number_of_inputs(model::GEBT)
    return 6*length(model.icol_pt) + 6*length(model.icol_beam)
end

number_of_parameters(model::GEBT) = 3*length(model.icol_pt) + 36*length(model.icol_beam)

inplaceness(::Type{<:GEBT}) = InPlace()

mass_matrix_type(::Type{<:GEBT}) = Linear()

state_jacobian_type(::Type{<:GEBT}) = Nonlinear()

input_jacobian_type(::Type{<:GEBT}) = Linear()

# --- Methods --- #

function get_rates!(dq, model::GEBT, q, r, p, t)

    # number of points and elements
    npoint = length(model.icol_pt)
    nelem = length(model.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, npoint, nelem, model.start, model.stop)

    # set prescribed and distributed displacements/loads from inputs
    prescribed, distributed = gxbeam_loads(r, npoint, nelem, model.isforce)

    # set origin, linear velocity, and angular velocity from inputs
    x0 = @SVector zeros(3)
    v0 = SVector(r[end-5], r[end-4], r[end-3])
    ω0 = SVector(r[end-2], r[end-1], r[end])

    # return mass matrix multiplied state rates
    return gxbeam_rates!(dq, q, assembly, prescribed, distributed,
        model.force_scaling, model.mass_scaling, model.irow_pt, model.irow_beam,
        model.irow_beam1, model.irow_beam2, model.icol_pt, model.icol_beam,
        x0, v0, ω0)
end

function get_mass_matrix!(M, model::GEBT, q, r, p, t)

    # number of points and elements
    npoint = length(model.icol_pt)
    nelem = length(model.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, npoint, nelem, model.start, model.stop)

    # extract scaling parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # extract model pointers
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

    # number of points and elements
    npoint = length(model.icol_pt)
    nelem = length(model.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, npoint, nelem, model.start, model.stop)

    # set prescribed and distributed displacements/loads from inputs
    prescribed, distributed = gxbeam_loads(r, npoint, nelem, model.isforce)

    # set origin, linear velocity, and angular velocity from inputs
    x0 = @SVector zeros(3)
    v0 = SVector(r[end-5], r[end-4], r[end-3])
    ω0 = SVector(r[end-2], r[end-1], r[end])

    # return jacobian of right hand side with respect to state variables
    return gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

# TODO: Update this function to include prescribed conditions and body velocities
# function get_input_jacobian(Jr, model::GEBT, q, r, p, t)
#
#     # extract system, assembly, and distributed loads
#     system = model.system
#     assembly = model.assembly
#     distributed = model.distributed
#
#     # extract start and stop of each beam element
#     elements = assembly.elements
#     start = assembly.start
#     stop = assembly.stop
#
#     # extract system constants and pointers
#     force_scaling = system.force_scaling
#     mass_scaling = system.mass_scaling
#     irow_pt = system.irow_pt
#
#     # get elements to which distributed loads are applied
#     element_indices = keys(distributed)
#
#     # return input jacobian linear map
#     return gxbeam_input_jacobian(q, r, elements, start, stop, force_scaling,
#         mass_scaling, irow_pt, element_indices)
# end

# --- Unit Testing Methods --- #

function get_lhs(model::GEBT, dq, q, r, p, t)

    # number of points and elements
    npoint = length(model.icol_pt)
    nelem = length(model.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, npoint, nelem, model.start, model.stop)

    # set prescribed and distributed displacements/loads from inputs
    prescribed, distributed = gxbeam_loads(r, npoint, nelem, model.isforce)

    # set origin, linear velocity, and angular velocity from inputs
    x0 = @SVector zeros(3)
    v0 = SVector(r[end-5], r[end-4], r[end-3])
    ω0 = SVector(r[end-2], r[end-1], r[end])

    return gxbeam_lhs(q, dq, assembly, prescribed, distributed, force_scaling,
        mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, x0, v0, ω0)
end

# --- Internal --- #

function gxbeam_assembly(p, np, ne, start, stop)
    points = [gxbeam_point(view(p, 3*(i-1) + 1 : 3*(i-1) + 3)) for i = 1:np]
    elements = [gxbeam_element(view(p, 3*np + 36*(i-1) + 1 : 3*np + 36*(i-1) + 36),
        start[i], stop[i]) for i = 1:ne)]
    return points, elements
end

gxbeam_point(p) = SVector(p[1], p[2], p[3])

function gxbeam_element(p, start, stop)
    # separate element parameters
    e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z, C11, C12, C13, C14, C15, C16,
        C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56,
        C66, μ, xm2, xm3, i22, i33, i23 = p
    # element length
    ΔL = points[stop[i]] - points[start[i]]
    # element location
    x = (points[start[i]] + points[stop[i]])/2
    # element compliance matrix
    C = @SMatrix [
        C11 C12 C13 C14 C15 C16;
        C12 C22 C23 C24 C25 C26;
        C13 C23 C33 C34 C35 C36;
        C14 C24 C34 C44 C45 C46;
        C15 C25 C35 C45 C55 C56;
        C16 C26 C36 C46 C56 C66
        ]
    # element inverse mass matrix
    minv = inv(@SMatrix [
            μ       0     0       0 μ*xm3 -μ*xm3;
            0       μ     0  -μ*xm3     0      0;
            0       0     μ   μ*xm2     0      0;
            0  -μ*xm3 μ*xm2 i22+i33     0      0;
         μ*xm3      0     0       0   i22   -i23;
        -μ*xm2      0     0       0  -i23    i33;
    ])
    # element triad
    Cab = @SMatrix [
        e1x e2x e3x;
        e1y e2y e3y;
        e1z e2z e3z;
    ]
    return Element(ΔL, x, C, minv, Cab)
end

function gxbeam_loads()

end

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
