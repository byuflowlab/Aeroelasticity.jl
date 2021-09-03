"""
    Initialize{GEBT} <: AbstractModel

Initial condition analysis model for the ``GEBT`` model.  This model may be used
to define a consistent set of initial conditions for use with the GEBT model.
State variables for this model are as defined by GXBeam.  Inputs correspond to
the external forces ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or
displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i},
\\theta_{z,i}`` applied to each node, followed by the distributed loads
``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}`` applied to each beam
element, followed by the linear and angular velocity of the system. Parameters
correspond to the location ``p_{x}, p_{y}, p_{z}`` of each node, followed by
each beam element's properties, followed by the initial displacement and velocity
of each beam element ``u_{x}, u_{y}, u_{z}, \dot{u}_x, \dot{u}_y, \dot{u}_z``.
Each beam element's properties are defined by a triad which defines the
orientation of the beam element ``e_{1,x}, e_{1,y}, e_{1,z}, e_{2,x}, e_{2,y},
e_{2,z}, e_{3,x}, e_{3,y}, e_{3,z}``, followed by the 21 independent entries of
the compliance matrix ``C_{11}, C_{12}, C_{13}, C_{14}, C_{15}, C_{16}, C_{22},
C_{23}, C_{24}, C_{25}, C_{26}, C_{33}, C_{34}, C_{35}, C_{36}, C_{44}, C_{45},
C_{46}, C_{55}, C_{56}, C_{66}``, followed by the beam element's inertial
properties ``\\mu, x_{m,2}, x_{m,3}, i_{22}, i_{33}, i_{23}``.
"""
struct Initialize{T}
    model::T
end

# --- Traits --- #

function number_of_states(model::Initialize)
    return number_of_states(model.model)
end

function number_of_inputs(model::Initialize)
    return number_of_inputs(model.model)
end

function number_of_parameters(model::Initialize{<:GEBT})
    return number_of_parameters(model.model) + 6*length(model.model.icol_elem)
end

inplaceness(::Type{Initialize{T}}) where T = inplaceness(T)

mass_matrix_type(::Type{Initialize{T}}) where T = mass_matrix_type(T)

state_jacobian_type(::Type{Initialize{T}}) where T = state_jacobian_type(T)

input_jacobian_type(::Type{Initialize{T}}) where T = input_jacobian_type(T)

# --- Methods --- #

function get_rates!(Mdq, model::Initialize{<:GEBT}, q, r, p, t)

    # number of points and elements
    np = length(model.model.icol_point)
    ne = length(model.model.icol_elem)

    # index of beginning and end of each beam element
    start = model.model.start
    stop = model.model.stop

    # scaling parameters
    force_scaling = model.model.force_scaling
    mass_scaling = model.model.mass_scaling

    # model pointers
    irow_point = model.model.irow_point
    irow_elem = model.model.irow_elem
    irow_elem1 = model.model.irow_elem1
    irow_elem2 = model.model.irow_elem2
    icol_point = model.model.icol_point
    icol_elem = model.model.icol_elem

    # location and degree of freedom of displacement constraints
    displacement = model.model.displacement

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # extract newly created beam elements
    elements = assembly.elements

    # construct point and element loads from inputs
    prescribed, distributed = gxbeam_loads(r, np, ne, displacement, elements)

    # origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector [r[end-5], r[end-4], r[end-3]]
    ω0 = @SVector [r[end-2], r[end-1], r[end]]

    # return mass matrix multiplied state rates
    return gxbeam_initialize_rates!(Mdq, q, assembly, prescribed, distributed,
        force_scaling, mass_scaling, irow_point, irow_elem,
        irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0)
end

function get_mass_matrix!(M, model::Initialize, q, r, p, t)

    # number of points and elements
    np = length(model.icol_point)
    ne = length(model.icol_elem)

    # index of beginning and end of each beam element
    start = model.start
    stop = model.stop

    # scaling parameters
    force_scaling = model.force_scaling
    mass_scaling = model.mass_scaling

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # extract model pointers
    irow_point = model.irow_point
    irow_elem = model.irow_elem
    irow_elem1 = model.irow_elem1
    irow_elem2 = model.irow_elem2
    icol_point = model.icol_point
    icol_elem = model.icol_elem

    # return mass matrix
    return gxbeam_mass_matrix!(M, q, assembly, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

# --- Performance Overloads --- #

function get_state_jacobian!(J, model::Initialize, q, r, p, t)

    # number of points and elements
    np = length(model.icol_point)
    ne = length(model.icol_elem)

    # index of beginning and end of each beam element
    start = model.start
    stop = model.stop

    # scaling parameters
    force_scaling = model.force_scaling
    mass_scaling = model.mass_scaling

    # model pointers
    irow_point = model.irow_point
    irow_elem = model.irow_elem
    irow_elem1 = model.irow_elem1
    irow_elem2 = model.irow_elem2
    icol_point = model.icol_point
    icol_elem = model.icol_elem

    # location and degree of freedom of displacement constraints
    displacement = model.displacement

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # extract newly created beam elements
    elements = assembly.elements

    # construct point and element loads from inputs
    prescribed, distributed = gxbeam_loads(r, np, ne, displacement, elements)

    # origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector [r[end-5], r[end-4], r[end-3]]
    ω0 = @SVector [r[end-2], r[end-1], r[end]]

    # return jacobian of right hand side with respect to state variables
    return gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_point, irow_elem,
        irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)
end

# TODO: Update this function to include prescribed conditions and body velocities
# function get_input_jacobian(Jr, model::Initialize, q, r, p, t)
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
#     irow_point = system.irow_point
#
#     # get elements to which distributed loads are applied
#     element_indices = keys(distributed)
#
#     # return input jacobian linear map
#     return gxbeam_input_jacobian(q, r, elements, start, stop, force_scaling,
#         mass_scaling, irow_point, element_indices)
# end

# --- Unit Testing Methods --- #

function get_lhs(model::Initialize, dq, q, r, p, t)

    # number of points and elements
    np = length(model.icol_point)
    ne = length(model.icol_elem)

    # index of beginning and end of each beam element
    start = model.start
    stop = model.stop

    # scaling parameters
    force_scaling = model.force_scaling
    mass_scaling = model.mass_scaling

    # model pointers
    irow_point = model.irow_point
    irow_elem = model.irow_elem
    irow_elem1 = model.irow_elem1
    irow_elem2 = model.irow_elem2
    icol_point = model.icol_point
    icol_elem = model.icol_elem

    # location and degree of freedom of displacement constraints
    displacement = model.displacement

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # extract newly created beam elements
    elements = assembly.elements

    # construct point and element loads from inputs
    prescribed, distributed = gxbeam_loads(r, np, ne, displacement, elements)

    # origin, linear velocity, and angular velocity
    x0 = @SVector zeros(3)
    v0 = @SVector [r[end-5], r[end-4], r[end-3]]
    ω0 = @SVector [r[end-2], r[end-1], r[end]]

    return gxbeam_lhs(q, dq, assembly, prescribed, distributed, force_scaling,
        mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point,
        icol_elem, x0, v0, ω0)
end

# --- Convenience Functions --- #

function set_states!(x, model::Initialize;
    u_e = nothing, theta_e = nothing, F_e = nothing, M_e = nothing,
    P_e = nothing, H_e = nothing, u_p = nothing, theta_p = nothing,
    F_p = nothing, M_p = nothing)

    icol_elem = model.icol_elem
    icol_point = model.icol_point

    displacement = model.displacement

    force_scaling = model.force_scaling
    mass_scaling = model.mass_scaling

    nelem = length(icol_elem)
    npoint = length(icol_point)

    for ielem = 1:nelem

        icol = icol_elem[ielem]

        if !isnothing(u_e)
            GXBeam.set_element_deflections!(x, icol, u_e[ielem])
        end

        if !isnothing(theta_e)
            GXBeam.set_element_rotations!(x, icol, theta_e[ielem])
        end

        if !isnothing(F_e)
            GXBeam.set_element_forces!(x, icol, F_e[ielem], force_scaling)
        end

        if !isnothing(M_e)
            GXBeam.set_element_moments!(x, icol, M_e[ielem], force_scaling)
        end

        if !isnothing(P_e)
            GXBeam.set_element_linear_momenta!(x, icol, P_e[ielem], mass_scaling)
        end

        if !isnothing(H_e)
            GXBeam.set_element_angular_momenta!(x, icol, H_e[ielem], mass_scaling)
        end
    end

    for ipoint = 1:npoint

        icol = icol_point[ipoint]

        prescribed_forces = SVector{6}(view(displacement, :, ipoint)) .== false

        if !isnothing(u_p)
            GXBeam.set_point_deflections!(x, icol, u_p[ipoint], prescribed_forces)
        end

        if !isnothing(theta_p)
            GXBeam.set_point_rotations!(x, icol, theta_p[ipoint], prescribed_forces)
        end

        if !isnothing(F_p)
            GXBeam.set_point_forces!(x, icol, F_p[ipoint], prescribed_forces, force_scaling)
        end

        if !isnothing(M_p)
            GXBeam.set_point_moments!(x, icol, M_p[ipoint], prescribed_forces, force_scaling)
        end
    end

    return x
end

function set_inputs!(y, model::Initialize; point_conditions, element_loads, V, Omega)

    np = length(model.icol_point)
    ne = length(model.icol_elem)

    for ip = 1:np
        y[6*(ip-1) + 1] = point_conditions[1,ip]
        y[6*(ip-1) + 2] = point_conditions[2,ip]
        y[6*(ip-1) + 3] = point_conditions[3,ip]
        y[6*(ip-1) + 4] = point_conditions[4,ip]
        y[6*(ip-1) + 5] = point_conditions[5,ip]
        y[6*(ip-1) + 6] = point_conditions[6,ip]
    end

    for ie = 1:ne
        y[6*np + 6*(ie-1) + 1] = element_loads[1,ie]
        y[6*np + 6*(ie-1) + 2] = element_loads[2,ie]
        y[6*np + 6*(ie-1) + 3] = element_loads[3,ie]
        y[6*np + 6*(ie-1) + 4] = element_loads[4,ie]
        y[6*np + 6*(ie-1) + 5] = element_loads[5,ie]
        y[6*np + 6*(ie-1) + 6] = element_loads[6,ie]
    end

    y[6*np + 6*ne + 1] = V[1]
    y[6*np + 6*ne + 2] = V[2]
    y[6*np + 6*ne + 3] = V[3]

    y[6*np + 6*ne + 4] = Omega[1]
    y[6*np + 6*ne + 5] = Omega[2]
    y[6*np + 6*ne + 6] = Omega[3]

    return y
end

function set_parameters!(p, model::Initialize; assembly)

    np = length(model.icol_point)
    ne = length(model.icol_elem)

    for ip = 1:np
        p[3*(ip-1) + 1] = assembly.points[ip][1]
        p[3*(ip-1) + 2] = assembly.points[ip][2]
        p[3*(ip-1) + 3] = assembly.points[ip][3]
    end

    for ie = 1:ne
        element = assembly.elements[ie]
        minv11, minv12, minv22 = element.minv11, element.minv12, element.minv22
        Msub = inv((@SMatrix [
            minv11[1,1] minv12[1,2] minv12[1,3];
            minv12[1,2] minv22[2,2] minv22[2,3];
            minv12[1,3] minv22[3,2] minv22[3,3];
            ]))
        p[3*np + 36*(ie-1) + 1] = element.Cab[1,1]
        p[3*np + 36*(ie-1) + 2] = element.Cab[2,1]
        p[3*np + 36*(ie-1) + 3] = element.Cab[3,1]
        p[3*np + 36*(ie-1) + 4] = element.Cab[1,2]
        p[3*np + 36*(ie-1) + 5] = element.Cab[2,2]
        p[3*np + 36*(ie-1) + 6] = element.Cab[3,2]
        p[3*np + 36*(ie-1) + 7] = element.Cab[1,3]
        p[3*np + 36*(ie-1) + 8] = element.Cab[2,3]
        p[3*np + 36*(ie-1) + 9] = element.Cab[3,3]
        p[3*np + 36*(ie-1) + 10] = element.C11[1,1]
        p[3*np + 36*(ie-1) + 11] = element.C11[1,2]
        p[3*np + 36*(ie-1) + 12] = element.C11[1,3]
        p[3*np + 36*(ie-1) + 13] = element.C12[1,1]
        p[3*np + 36*(ie-1) + 14] = element.C12[1,2]
        p[3*np + 36*(ie-1) + 15] = element.C12[1,3]
        p[3*np + 36*(ie-1) + 16] = element.C11[2,2]
        p[3*np + 36*(ie-1) + 17] = element.C11[2,3]
        p[3*np + 36*(ie-1) + 18] = element.C12[2,1]
        p[3*np + 36*(ie-1) + 19] = element.C12[2,2]
        p[3*np + 36*(ie-1) + 20] = element.C12[2,3]
        p[3*np + 36*(ie-1) + 21] = element.C11[3,3]
        p[3*np + 36*(ie-1) + 22] = element.C12[3,1]
        p[3*np + 36*(ie-1) + 23] = element.C12[3,2]
        p[3*np + 36*(ie-1) + 24] = element.C12[3,3]
        p[3*np + 36*(ie-1) + 25] = element.C22[1,1]
        p[3*np + 36*(ie-1) + 26] = element.C22[1,2]
        p[3*np + 36*(ie-1) + 27] = element.C22[1,3]
        p[3*np + 36*(ie-1) + 28] = element.C22[2,2]
        p[3*np + 36*(ie-1) + 29] = element.C22[2,3]
        p[3*np + 36*(ie-1) + 30] = element.C22[3,3]
        p[3*np + 36*(ie-1) + 31] =  Msub[1,1] # μ
        p[3*np + 36*(ie-1) + 32] = -Msub[1,3]/Msub[1,1] # xm2
        p[3*np + 36*(ie-1) + 33] =  Msub[1,2]/Msub[1,1] # xm3
        p[3*np + 36*(ie-1) + 34] =  Msub[2,2] # i22
        p[3*np + 36*(ie-1) + 35] =  Msub[3,3] # i33
        p[3*np + 36*(ie-1) + 36] = -Msub[2,3] # i23
    end

    return p
end

function separate_states(model::Initialize, x;
    point_conditions = fill(NaN, 6, length(model.icol_point)))

    TF = eltype(x)
    np = length(model.icol_point)
    ne = length(model.icol_elem)

    u_p = Vector{SVector{3,TF}}(undef, np)
    theta_p = Vector{SVector{3,TF}}(undef, np)
    F_p = Vector{SVector{3,TF}}(undef, np)
    M_p = Vector{SVector{3,TF}}(undef, np)

    u_e = Vector{SVector{3,TF}}(undef, ne)
    theta_e = Vector{SVector{3,TF}}(undef, ne)
    F_e = Vector{SVector{3,TF}}(undef, ne)
    M_e = Vector{SVector{3,TF}}(undef, ne)
    P_e = Vector{SVector{3,TF}}(undef, ne)
    H_e = Vector{SVector{3,TF}}(undef, ne)

    for ip = 1:np

        icol = model.icol_point[ip]

        prescribed_forces = SVector{6}(view(model.displacement, :, ip)) .== false

        # get the displacement and rotations of the point
        u_p[ip] = SVector(ifelse(prescribed_forces[1], x[icol  ], point_conditions[1,ip]),
                    ifelse(prescribed_forces[2], x[icol+1], point_conditions[2,ip]),
                    ifelse(prescribed_forces[3], x[icol+2], point_conditions[3,ip]))
        theta_p[ip] = SVector(ifelse(prescribed_forces[4], x[icol+3], point_conditions[4,ip]),
                    ifelse(prescribed_forces[5], x[icol+4], point_conditions[5,ip]),
                    ifelse(prescribed_forces[6], x[icol+5], point_conditions[6,ip]))

        # overwrite external forces/moments with solved for forces/moments
        F_p[ip] = SVector(ifelse(prescribed_forces[1], point_conditions[1,ip], x[icol  ] * model.force_scaling),
                    ifelse(prescribed_forces[2], point_conditions[2,ip], x[icol+1] * model.force_scaling),
                    ifelse(prescribed_forces[3], point_conditions[3,ip], x[icol+2] * model.force_scaling))
        M_p[ip] = SVector(ifelse(prescribed_forces[4], point_conditions[4,ip], x[icol+3] * model.force_scaling),
                    ifelse(prescribed_forces[5], point_conditions[5,ip], x[icol+4] * model.force_scaling),
                    ifelse(prescribed_forces[6], point_conditions[6,ip], x[icol+5] * model.force_scaling))

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(theta_p[ip])
        theta_p[ip] *= scaling

    end

    for ie = 1:ne

        icol = model.icol_elem[ie]

        u_e[ie] = SVector(x[icol], x[icol+1], x[icol+2])
        theta_e[ie] = SVector(x[icol+3], x[icol+4], x[icol+5])
        F_e[ie] = SVector(x[icol+6], x[icol+7], x[icol+8]) .* model.force_scaling
        M_e[ie] = SVector(x[icol+9], x[icol+10], x[icol+11]) .* model.force_scaling
        P_e[ie] = SVector(x[icol+12], x[icol+13], x[icol+14]) .* model.mass_scaling
        H_e[ie] = SVector(x[icol+15], x[icol+16], x[icol+17]) .* model.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(theta_e[ie])
        theta_e[ie] *= scaling

    end

    return (u_p = u_p, theta_p = theta_p, F_p = F_p, M_p = M_p,
        u_e = u_e, theta_e = theta_e, F_e = F_e, M_e = M_e, P_e = P_e, H_e = H_e)
end

function separate_inputs(model::Initialize, y)

    np = length(model.icol_point)
    ne = length(model.icol_elem)

    point_conditions = reshape(view(y, 1 : 6*np), 6, np)

    element_loads = reshape(view(y, 6*np + 1 : 6*np + 6*ne), 6, ne)

    V = view(y, 6*np + 6*ne + 1 : 6*np + 6*ne + 3)

    Omega = view(y, 6*np + 6*ne + 4 : 6*np + 6*ne + 6)

    return (point_conditions = point_conditions, element_loads = element_loads,
        V = V, Omega = Omega)
end

function separate_parameters(model::Initialize, p)
    np = length(model.icol_point)
    ne = length(model.icol_elem)
    start = model.start
    stop = model.stop
    return (assembly = gxbeam_assembly(p, np, ne, start, stop),)
end

# --- Internal --- #

function gxbeam_assembly(p, np, ne, start, stop)
    points = [gxbeam_point(view(p, 3*(ip-1) + 1 : 3*(ip-1) + 3)) for ip = 1:np]
    elements = [gxbeam_element(view(p, 3*np + 36*(ie-1) + 1 : 3*np + 36*(ie-1) + 36),
        points, start[ie], stop[ie]) for ie = 1:ne]
    return GXBeam.Assembly(points, start, stop, elements)
end

gxbeam_point(p) = SVector(p[1], p[2], p[3])

function gxbeam_element(p, points, start, stop)
    # separate element parameters
    e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z, C11, C12, C13, C14, C15, C16,
        C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56,
        C66, μ, xm2, xm3, i22, i33, i23 = p
    # element length
    ΔL = norm(points[stop] - points[start])
    # element location
    x = (points[start] + points[stop])/2
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
            μ       0     0       0 μ*xm3 -μ*xm2;
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

    return GXBeam.Element(ΔL, x, C, minv, Cab)
end

function gxbeam_loads(r, np, ne, d, elements)
    prescribed = Dict(ip => gxbeam_point_load(view(r,
        6*(ip-1) + 1 : 6*(ip-1) + 6), view(d, :, ip)) for ip = 1:np)
    distributed = Dict(ie => gxbeam_distributed_load(view(r,
        6*np + 6*(ie-1) + 1 : 6*np + 6*(ie-1) + 6), elements[ie].L) for ie = 1:ne)
    return prescribed, distributed
end

function gxbeam_point_load(r, d)
    displacement = SVector{6}(d)
    force = displacement .== false
    value = SVector{6}(r)
    follower = @SVector zeros(6)
    return PrescribedConditions(force, value, follower)
end

function gxbeam_distributed_load(r, ΔL)
    f = SVector(r[1], r[2], r[3])
    m = SVector(r[4], r[5], r[6])
    f1 = f2 = ΔL*f/2
    m1 = m2 = ΔL*m/2
    f1_follower = f2_follower = @SVector zeros(3)
    m1_follower = m2_follower = @SVector zeros(3)
    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower,
        m2_follower)
end

function gxbeam_lhs(q, dq, assembly, prescribed, distributed,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
    icol_point, icol_elem, x0, v0, ω0)

    steady_residual = similar(dq)
    dynamic_residual = similar(dq)

    GXBeam.steady_state_system_residual!(steady_residual, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1,
        irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    GXBeam.dynamic_system_residual!(dynamic_residual, q, dq, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1,
        irow_elem2, icol_point, icol_elem, x0, v0, ω0)



    j! = (J, x) -> system_jacobian!(J, x, assembly, pcond, dload, force_scaling,
        mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, u0, theta0, udot0, thetadot0)

    return steady_residual - dynamic_residual
end

function gxbeam_rhs!(out, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_point, irow_elem,
    irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    GXBeam.system_residual!(F, x, assembly, pcond, dload, force_scaling,
        mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, u0, theta0, udot0, thetadot0)

    return out
end

function gxbeam_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_point,
    irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    M .= 0

    # mass matrix is GXBeam mass matrix, moved to LHS
    GXBeam.system_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_point,
        irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    M .*= -1

    return M
end

function gxbeam_rates!(Mdu, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_point, irow_elem,
    irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    gxbeam_rhs!(Mdu, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_point, irow_elem,
        irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    return Mdu
end

function gxbeam_state_jacobian!(K, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_point, irow_elem,
    irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    K .= 0

    # jacobian of mass matrix multiplied rates is equal to steady state jacobian
    GXBeam.steady_state_system_jacobian!(K, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1,
        irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    return K
end

# function gxbeam_input_jacobian(q, r, elements, start, stop, force_scaling, mass_scaling,
#     irow_point, element_indices)
#
#     f! = (y, x) -> gxbeam_jac_vec!(y, x, r, elements, start, stop, force_scaling,
#         mass_scaling, irow_point, element_indices)
#     M = length(q)
#     N = length(r)
#
#     return LinearMap(f!, M, N; ismutating = true)
# end
#
# # input jacobian vector product
# function gxbeam_jac_vec!(y, x, r, elements, start, stop, force_scaling, mass_scaling,
#     irow_point, element_indices)
#
#     nelem = length(start)
#
#     # initialize output array
#     y .= 0
#
#     # current input array index
#     ir = 0
#
#     # distributed loads
#     for ielem = 1:nelem
#         if ielem in element_indices
#             # get point indices
#             irow_p1 = irow_point[start[ielem]]
#             irow_p2 = irow_point[stop[ielem]]
#             # initialize element load jacobians
#             f_r = Diagonal((@SVector ones(3)))
#             m_r = Diagonal((@SVector ones(3)))
#             # jacobians after integration
#             # (assume distributed loads are constant on each element)
#             f1_r = f2_r = elements[ielem].L*f_r/2
#             m1_r = m2_r = elements[ielem].L*m_r/2
#             # calculate jacobian vector product
#             xf = SVector(x[ir+1], x[ir+2], x[ir+3])
#             xm = SVector(x[ir+4], x[ir+5], x[ir+6])
#             y[irow_p1:irow_p1+2] .+= -f1_r * xf ./ force_scaling
#             y[irow_p1+3:irow_p1+5] .+= -m1_r * xm ./ force_scaling
#             y[irow_p2:irow_p2+2] .+= -f2_r * xf ./ force_scaling
#             y[irow_p2+3:irow_p2+5] .+= -m2_r * xm ./ force_scaling
#             # move to next set of distributed loads
#             ir += 6
#         end
#     end
#
#     return y
# end
#
# # update distributed loads using inputs
# function update_distributed(distributed, elements, r::AbstractVector{T}) where T
#     new_distributed = Dict{Int, DistributedLoads{T}}()
#     for (key, value) in distributed
#         new_distributed[key] = DistributedLoads{T}(value)
#     end
#     update_distributed!(new_distributed, elements, r)
#     return new_distributed
# end
#
# # update distributed loads using inputs
# function update_distributed!(distributed, elements, r)
#     distributed_points = keys(distributed)
#     largest_point = maximum(distributed_points)
#     ir = 0
#     for ielem in 1:largest_point
#         if ielem in distributed_points
#             # extract loads on this element
#             f = SVector(r[ir+1], r[ir+2], r[ir+3])
#             m = SVector(r[ir+4], r[ir+5], r[ir+6])
#             # get integrated distributed loads
#             f1 = f2 = elements[ielem].L*f/2
#             m1 = m2 = elements[ielem].L*m/2
#             # keep the same follower loads
#             f1_follower = distributed[ielem].f1_follower
#             f2_follower = distributed[ielem].f2_follower
#             m1_follower = distributed[ielem].m1_follower
#             m2_follower = distributed[ielem].m2_follower
#             # update distributed loads
#             distributed[ielem] = DistributedLoads(f1, f2, m1, m2, f1_follower,
#                 f2_follower, m1_follower, m2_follower)
#             ir += 6
#         end
#     end
#
#     return distributed
# end
