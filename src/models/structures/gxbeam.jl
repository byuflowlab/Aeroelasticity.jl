"""
    GXBeamAssembly

Model which describes the behavior of an assembly of beam elements, as modeled by the 
`GXBeam` package.  State variables are as defined by GXBeam.  Inputs correspond to the 
external forces ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or
displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i},
\\theta_{z,i}`` applied to each node, followed by the distributed loads
``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}`` applied to each beam
element, followed by the point mass properties ``m, p, I_{11}, I_{22}, I_{33}, I_{12}, 
I_{13}, I_{23}`` corresponding to each beam element, followed by the gravity vector, linear 
velocity, angular velocity, linear acceleration, and angular acceleration of the system. 
Parameters correspond to the location ``p_{x}, p_{y}, p_{z}`` of each node followed by each 
beam element's properties. Each beam element's properties are defined by a triad which defines 
the orientation of the beam element ``e_{1,x}, e_{1,y}, e_{1,z}, e_{2,x}, e_{2,y}, e_{2,z}, 
e_{3,x}, e_{3,y}, e_{3,z}``, followed by the 21 independent entries of the compliance matrix 
``C_{11}, C_{12}, C_{13}, C_{14}, C_{15}, C_{16}, C_{22}, C_{23}, C_{24}, C_{25}, C_{26}, 
C_{33}, C_{34}, C_{35}, C_{36}, C_{44}, C_{45}, C_{46}, C_{55}, C_{56}, C_{66}``, followed 
by the beam element's inertial properties ``\\mu, x_{m,2}, x_{m,3}, i_{22}, i_{33}, i_{23}``
"""
struct GXBeamAssembly{TF,TI,TC,TB}
    # beam element connectivity
    start::TC
    stop::TC
    # location and degree of freedom of displacement constraints
    displacement::TB
    # scaling parameters for the system
    force_scaling::TF
    # indices for accessing governing equations for each point and beam element
    irow_point::TI
    irow_elem::TI
    irow_elem1::TI
    irow_elem2::TI
    # indices for accessing state variables for each point and beam element
    icol_point::TI
    icol_elem::TI
end

# --- Constructors --- #

"""
    GXBeamAssembly(start, stop, displacement; kwargs...)

Construct a geometrically exact beam theory structural model with beam elements
which extend from the point indices in `start` to the point indices in `stop` and 
points with prescribed displacements as specified in `displacement`.

# Arguments
 - `start`: Vector containing point index where each beam element starts
 - `stop`: Vector containing point index where each beam element stops
 - `displacement`: Boolean matrix indicating the point index and degree of
    freedom of the system's displacement constraints.  Rows correspond to the
    degrees of freedom ``x, y, z, \\theta_x, \\theta_y, \\theta_z`` and
    columns correspond to the point indices ``1, 2, \\dots, N_p`` where ``N_p``
    is the total number of points.

# Keyword Arguments
 - `force_scaling = 1.0`: Factor used to scale forces and moments internally
"""
function GXBeamAssembly(start, stop, displacement; force_scaling=1.0)

    static = false

    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem =
        GXBeam.system_indices(start, stop, static)

    return GXBeamAssembly(start, stop, displacement, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

"""
    GXBeamAssembly(assembly, prescribed_conditions; kwargs...)

Construct a geometrically exact beam theory structural model with connectivity
as specified in `assembly` and displacement constraints as specified in 
`prescribed_conditions`.

# Keyword Arguments
 - `force_scaling = 1.0`: Factor used to scale forces and moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    compliance matrix.
"""
function GXBeamAssembly(assembly, prescribed_conditions;
    force_scaling = GXBeam.default_force_scaling(assembly))

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    start = assembly.start
    stop = assembly.stop

    static = false

    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem =
        GXBeam.system_indices(start, stop, static)

    displacement = zeros(Bool, 6, npoint)
    for key in keys(prescribed_conditions)
        displacement[:,key] .= prescribed_conditions[key].isforce .== false
    end

    return GXBeamAssembly(start, stop, displacement, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

# --- Submodel Creation --- #

function Submodel(model::GXBeamAssembly)
    
    # unpack constants
    start = model.start
    stop = model.stop
    displacement = model.displacement
    force_scaling = model.force_scaling
    irow_point = model.irow_point
    irow_elem = model.irow_elem
    irow_elem1 = model.irow_elem1
    irow_elem2 = model.irow_elem2
    icol_point = model.icol_point
    icol_elem = model.icol_elem

    # number of points and elements
    npoint = length(icol_point)
    nelem = length(icol_elem)

    # residual function
    fresid = (resid, dx, x, y, p, t) -> gxbeam_residual!(resid, dx, x, y, p, t; 
        start, stop, displacement, force_scaling, irow_point, irow_elem, 
        irow_elem1, irow_elem2, icol_point, icol_elem)

    # number of states, inputs, and parameters
    nx = 6*npoint + 18*nelem
    ny = 6*npoint + 6*nelem + 10*nelem + 15
    np = 3*npoint + 36*nelem

    # jacobian definitions
    ratejac = Linear((J, dx, x, y, p, t) -> gxbeam_rate_jacobian!(J, dx, x, y, p, t; 
        start, stop, displacement, force_scaling, irow_point, irow_elem, irow_elem1, 
        irow_elem2, icol_point, icol_elem))
    statejac = Nonlinear((J, dx, x, y, p, t) -> gxbeam_state_jacobian!(J, dx, x, y, p, t; 
        start, stop, displacement, force_scaling, irow_point, irow_elem, irow_elem1, 
        irow_elem2, icol_point, icol_elem))
    inputjac = Nonlinear() # TODO: define inputjac function
    paramjac = Nonlinear() # TODO: define paramjac function
    tgrad = Zeros()

    # convenience functions for setting states, inputs, and parameters
    setstate = (x; kwargs...) -> gxbeam_set_states!(x, displacement, icol_point, icol_elem, 
        force_scaling; kwargs...)

    setinput = (y; kwargs...) -> gxbeam_set_inputs!(y, displacement, icol_point, icol_elem; 
        kwargs...) 

    setparam = (p; kwargs...) -> gxbeam_set_parameters!(p, icol_point, icol_elem; kwargs...)

    # convenience functions for separating states, inputs, and parameters
    sepstate = (x) -> gxbeam_separate_states(x, displacement, icol_point, icol_elem, force_scaling)
    sepinput = (y) -> gxbeam_separate_inputs(y, displacement, icol_point, icol_elem)
    sepparam = (p) -> gxbeam_separate_parameters(p, start, stop, icol_point, icol_elem)

    # model definition
    return Submodel{true}(fresid, nx, ny, np;
        ratejac = ratejac,
        statejac = statejac,
        inputjac = inputjac,
        paramjac = paramjac,
        tgrad = tgrad,
        setstate = setstate,
        setinput = setinput,
        setparam = setparam,
        sepstate = sepstate,
        sepinput = sepinput,
        sepparam = sepparam,
    )
end

# --- Internal Methods for this Model --- #

# residual function
function gxbeam_residual!(resid, dx, x, y, p, t; start, stop, displacement, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    # number of points and elements
    np = length(icol_point)
    ne = length(icol_elem)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # construct point loads, distributed loads, and point masses
    prescribed, distributed, point_masses = gxbeam_loads(y, np, ne, displacement, assembly.elements)

    # gravity vector
    gvec = SVector(y[end-14], y[end-13], y[end-12])

    # linear/angular velocity and acceleration
    x0 = @SVector zeros(3)
    v0 = SVector(y[end-11], y[end-10], y[end-9])
    ω0 = SVector(y[end-8], y[end-7], y[end-6])
    a0 = SVector(y[end-5], y[end-4], y[end-3])
    α0 = SVector(y[end-2], y[end-1], y[end])

    # return mass matrix multiplied state rates
    return GXBeam.dynamic_system_residual!(resid, dx, x, assembly, prescribed, distributed, 
        point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, 
        icol_point, icol_elem, x0, v0, ω0, a0, α0)
end

# rate jacobian function
function gxbeam_rate_jacobian!(J, dx, x, y, p, t; start, stop, displacement, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    # zero out jacobian entries
    J .= 0

    # number of points and elements
    np = length(icol_point)
    ne = length(icol_elem)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # construct point loads, distributed loads, and point masses
    prescribed, distributed, point_masses = gxbeam_loads(y, np, ne, displacement, assembly.elements)

    # populate mass matrix entries
    return GXBeam.system_mass_matrix!(J, x, assembly, point_masses, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

# state jacobian function
function gxbeam_state_jacobian!(J, dx, x, y, p, t; start, stop, displacement, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    # zero out jacobian entries
    J .= 0

    # number of points and elements
    np = length(icol_point)
    ne = length(icol_elem)

    # construct assembly from parameters
    assembly = gxbeam_assembly(p, np, ne, start, stop)

    # construct point and element loads from inputs
    prescribed, distributed, point_masses = gxbeam_loads(y, np, ne, displacement, assembly.elements)

    # gravity vector
    gvec = SVector(y[end-14], y[end-13], y[end-12])

    # linear/angular velocity and acceleration
    x0 = @SVector zeros(3)
    v0 = SVector(y[end-11], y[end-10], y[end-9])
    ω0 = SVector(y[end-8], y[end-7], y[end-6])
    a0 = SVector(y[end-5], y[end-4], y[end-3])
    α0 = SVector(y[end-2], y[end-1], y[end])

    # return jacobian of right hand side with respect to state variables
    return GXBeam.dynamic_system_jacobian!(J, dx, x, assembly, prescribed, distributed, 
        point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, 
        icol_point, icol_elem, x0, v0, ω0, a0, α0)
end

# function for creating an assembly from the provided parameters
function gxbeam_assembly(p, np, ne, start, stop)
    points = [gxbeam_point(view(p, 3*(ip-1) + 1 : 3*(ip-1) + 3)) for ip = 1:np]
    elements = [gxbeam_element(view(p, 3*np + 36*(ie-1) + 1 : 3*np + 36*(ie-1) + 36),
        points, start[ie], stop[ie]) for ie = 1:ne]
    return GXBeam.Assembly(points, start, stop, elements)
end

# function for creating a point
gxbeam_point(p) = SVector(p[1], p[2], p[3])

# function for creating a beam element
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
    # element mass matrix
    mass = @SMatrix [
            μ       0     0       0 μ*xm3 -μ*xm2;
            0       μ     0  -μ*xm3     0      0;
            0       0     μ   μ*xm2     0      0;
            0  -μ*xm3 μ*xm2 i22+i33     0      0;
         μ*xm3      0     0       0   i22   -i23;
        -μ*xm2      0     0       0  -i23    i33;
    ]
    # element triad
    Cab = @SMatrix [
        e1x e2x e3x;
        e1y e2y e3y;
        e1z e2z e3z;
    ]

    return GXBeam.Element(ΔL, x, C, mass, Cab)
end

# define prescribed conditions, distributed loads, and point masses
function gxbeam_loads(y, np, ne, d, elements)
    prescribed = Dict(ip => gxbeam_prescribed_condition(view(y,
        6*(ip-1) + 1 : 6*(ip-1) + 6), view(d, :, ip)) for ip = 1:np)
    distributed = Dict(ie => gxbeam_distributed_load(view(y,
        6*np + 6*(ie-1) + 1 : 6*np + 6*(ie-1) + 6), elements[ie].L) for ie = 1:ne)
    point_masses = Dict(ie => gxbeam_point_mass(view(y,
        6*np + 6*ne + 10*(ie-1) + 1 : 6*np + 6*ne + 10*(ie-1) + 10)) for ie = 1:ne)    
    return prescribed, distributed, point_masses
end

# define prescribed conditions at a point
function gxbeam_prescribed_condition(y, d)
    displacement = SVector{6}(d)
    force = displacement .== false
    value = SVector{6}(y)
    follower = @SVector zeros(6)
    return PrescribedConditions(force, value, follower)
end

# define distributed loads on an element
function gxbeam_distributed_load(y, ΔL)
    # prescribed distributed loads
    f = SVector(y[1], y[2], y[3])
    m = SVector(y[4], y[5], y[6])
    # non-follower distributed loads
    f1 = f2 = ΔL*f/2
    m1 = m2 = ΔL*m/2
    # follower distributed loads
    f1_follower = f2_follower = @SVector zeros(3)
    m1_follower = m2_follower = @SVector zeros(3)
    # distributed loads
    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower,
        m2_follower)
end

# define point mass attached to an element
function gxbeam_point_mass(y)
    # point mass
    m = y[1]
    # point location
    p = SVector(y[2], y[3], y[4])
    # point inertia
    I11 = y[5]
    I22 = y[6]
    I33 = y[7]
    I12 = y[8]
    I13 = y[9]
    I23 = y[10]
    I = @SMatrix [I11 I12 I13; I12 I22 I23; I13 I23 I33]
    # point mass
    return PointMass(m, p, I)
end

# convenience function for defining this model's state vector
function gxbeam_set_states!(x, displacement, icol_point, icol_elem, force_scaling;
    u_e = nothing, theta_e = nothing, F_e = nothing, M_e = nothing,
    V_e = nothing, Omega_e = nothing, u_p = nothing, theta_p = nothing,
    F_p = nothing, M_p = nothing)

    nelem = length(icol_elem)
    npoint = length(icol_point)

    for ielem = 1:nelem
        icol = icol_elem[ielem]
        if !isnothing(u_e)
            GXBeam.set_element_deflection!(x, icol, u_e[ielem])
        end
        if !isnothing(theta_e)
            GXBeam.set_element_rotation!(x, icol, theta_e[ielem])
        end
        if !isnothing(F_e)
            GXBeam.set_element_forces!(x, icol, F_e[ielem], force_scaling)
        end
        if !isnothing(M_e)
            GXBeam.set_element_moments!(x, icol, M_e[ielem], force_scaling)
        end
        if !isnothing(V_e)
            GXBeam.set_element_linear_velocity!(x, icol, V_e[ielem])
        end
        if !isnothing(Omega_e)
            GXBeam.set_element_angular_velocity!(x, icol, Omega_e[ielem])
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

# convenience function for defining this model's input vector
function gxbeam_set_inputs!(y, displacement, icol_point, icol_elem; 
    u_p = nothing, theta_p = nothing, F_p = nothing, M_p = nothing, 
    f_d = nothing, m_d = nothing, 
    m_m = nothing, p_m = nothing, I_m = nothing,
    gvec = nothing, v = nothing, omega = nothing, a = nothing, alpha = nothing)

    np = length(icol_point)
    ne = length(icol_elem)

    for ip = 1:np
        prescribed_forces = SVector{6}(view(displacement, :, ip)) .== false
        for i = 1:3
            if prescribed_forces[i]
                if !isnothing(F_p)
                    y[6*(ip-1) + i] = F_p[ip][i]
                end
            else
                if !isnothing(u_p)
                    y[6*(ip-1) + i] = u_p[ip][i]
                end
            end
        end
        for i = 1:3
            if prescribed_forces[i]
                if !isnothing(M_p)
                    y[6*(ip-1) + 3 + i] = M_p[ip][i]
                end
            else
                if !isnothing(theta_p)
                    y[6*(ip-1) + 3 + i] = theta_p[ip][i]
                end
            end
        end
    end

    for ie = 1:ne
        if !isnothing(f_d)
            for i = 1:3
                y[6*np + 6*(ie-1) + i] = f_d[ie][i]
            end
        end
        if !isnothing(m_d)
            for i = 1:3
                y[6*np + 6*(ie-1) + 3 + i] = m_d[ie][i]
            end
        end
    end

    for ie = 1:ne
        if !isnothing(m_m)
            y[6*np + 6*ne + 10*(ie-1) + 1] = m_m[ie]
        end
        if !isnothing(p_m)
            y[6*np + 6*ne + 10*(ie-1) + 2] = p_m[ie][1]
            y[6*np + 6*ne + 10*(ie-1) + 3] = p_m[ie][2]
            y[6*np + 6*ne + 10*(ie-1) + 4] = p_m[ie][3]
        end
        if !isnothing(I_m)
            y[6*np + 6*ne + 10*(ie-1) + 5] = I_m[ie][1,1]
            y[6*np + 6*ne + 10*(ie-1) + 6] = I_m[ie][2,2]
            y[6*np + 6*ne + 10*(ie-1) + 7] = I_m[ie][3,3]
            y[6*np + 6*ne + 10*(ie-1) + 8] = I_m[ie][1,2]
            y[6*np + 6*ne + 10*(ie-1) + 9] = I_m[ie][1,3]
            y[6*np + 6*ne + 10*(ie-1) + 10] = I_m[ie][2,3]
        end
    end

    if !isnothing(gvec)
        y[6*np + 6*ne + 10*ne + 1] = gvec[1]
        y[6*np + 6*ne + 10*ne + 2] = gvec[2]
        y[6*np + 6*ne + 10*ne + 3] = gvec[3]
    end

    if !isnothing(v)
        y[6*np + 6*ne + 10*ne + 4] = v[1]
        y[6*np + 6*ne + 10*ne + 5] = v[2]
        y[6*np + 6*ne + 10*ne + 6] = v[3]
    end

    if !isnothing(omega)
        y[6*np + 6*ne + 10*ne + 7] = omega[1]
        y[6*np + 6*ne + 10*ne + 8] = omega[2]
        y[6*np + 6*ne + 10*ne + 9] = omega[3]
    end

    if !isnothing(a)
        y[6*np + 6*ne + 10*ne + 10] = a[1]
        y[6*np + 6*ne + 10*ne + 11] = a[2]
        y[6*np + 6*ne + 10*ne + 12] = a[3]
    end

    if !isnothing(alpha)
        y[6*np + 6*ne + 10*ne + 13] = alpha[1]
        y[6*np + 6*ne + 10*ne + 14] = alpha[2]
        y[6*np + 6*ne + 10*ne + 15] = alpha[3]
    end

    return y
end

# convenience function for defining this model's parameter vector
function gxbeam_set_parameters!(p, icol_point, icol_elem; assembly)

    np = length(icol_point)
    ne = length(icol_elem)

    for ip = 1:np
        p[3*(ip-1) + 1] = assembly.points[ip][1]
        p[3*(ip-1) + 2] = assembly.points[ip][2]
        p[3*(ip-1) + 3] = assembly.points[ip][3]
    end

    for ie = 1:ne
        element = assembly.elements[ie]
        p[3*np + 36*(ie-1) + 1] = element.Cab[1,1]
        p[3*np + 36*(ie-1) + 2] = element.Cab[2,1]
        p[3*np + 36*(ie-1) + 3] = element.Cab[3,1]
        p[3*np + 36*(ie-1) + 4] = element.Cab[1,2]
        p[3*np + 36*(ie-1) + 5] = element.Cab[2,2]
        p[3*np + 36*(ie-1) + 6] = element.Cab[3,2]
        p[3*np + 36*(ie-1) + 7] = element.Cab[1,3]
        p[3*np + 36*(ie-1) + 8] = element.Cab[2,3]
        p[3*np + 36*(ie-1) + 9] = element.Cab[3,3]
        p[3*np + 36*(ie-1) + 10] = element.compliance[1,1]
        p[3*np + 36*(ie-1) + 11] = element.compliance[1,2]
        p[3*np + 36*(ie-1) + 12] = element.compliance[1,3]
        p[3*np + 36*(ie-1) + 13] = element.compliance[1,4]
        p[3*np + 36*(ie-1) + 14] = element.compliance[1,5]
        p[3*np + 36*(ie-1) + 15] = element.compliance[1,6]
        p[3*np + 36*(ie-1) + 16] = element.compliance[2,2]
        p[3*np + 36*(ie-1) + 17] = element.compliance[2,3]
        p[3*np + 36*(ie-1) + 18] = element.compliance[2,4]
        p[3*np + 36*(ie-1) + 19] = element.compliance[2,5]
        p[3*np + 36*(ie-1) + 20] = element.compliance[2,6]
        p[3*np + 36*(ie-1) + 21] = element.compliance[3,3]
        p[3*np + 36*(ie-1) + 22] = element.compliance[3,4]
        p[3*np + 36*(ie-1) + 23] = element.compliance[3,5]
        p[3*np + 36*(ie-1) + 24] = element.compliance[3,6]
        p[3*np + 36*(ie-1) + 25] = element.compliance[4,4]
        p[3*np + 36*(ie-1) + 26] = element.compliance[4,5]
        p[3*np + 36*(ie-1) + 27] = element.compliance[4,6]
        p[3*np + 36*(ie-1) + 28] = element.compliance[5,5]
        p[3*np + 36*(ie-1) + 29] = element.compliance[5,6]
        p[3*np + 36*(ie-1) + 30] = element.compliance[6,6]
        p[3*np + 36*(ie-1) + 31] =  element.mass[1,1] # μ
        p[3*np + 36*(ie-1) + 32] =  element.mass[3,4]/element.mass[1,1] # xm2
        p[3*np + 36*(ie-1) + 33] =  element.mass[1,5]/element.mass[1,1] # xm3
        p[3*np + 36*(ie-1) + 34] =  element.mass[5,5] # i22
        p[3*np + 36*(ie-1) + 35] =  element.mass[6,6] # i33
        p[3*np + 36*(ie-1) + 36] = -element.mass[5,6] # i23
    end

    return p
end

# convenience function for separating this model's state vector
function gxbeam_separate_states(x, displacement, icol_point, icol_elem, force_scaling)
  
    # floating point type
    TF = eltype(x)

    # number of points and elements
    np = length(icol_point)
    ne = length(icol_elem)

    # separate point state variables
    u_p = Vector{SVector{3,TF}}(undef, np)
    theta_p = Vector{SVector{3,TF}}(undef, np)
    F_p = Vector{SVector{3,TF}}(undef, np)
    M_p = Vector{SVector{3,TF}}(undef, np)

    for ip = 1:np

        icol = icol_point[ip]

        prescribed_forces = SVector{6}(view(displacement, :, ip)) .== false

        # get the displacement and rotations of the point
        u_p[ip] = SVector(ifelse(prescribed_forces[1], x[icol  ], NaN),
                    ifelse(prescribed_forces[2], x[icol+1], NaN),
                    ifelse(prescribed_forces[3], x[icol+2], NaN))
        theta_p[ip] = SVector(ifelse(prescribed_forces[4], x[icol+3], NaN),
                    ifelse(prescribed_forces[5], x[icol+4], NaN),
                    ifelse(prescribed_forces[6], x[icol+5], NaN))

        # overwrite external forces/moments with solved for forces/moments
        F_p[ip] = SVector(ifelse(prescribed_forces[1], NaN, x[icol  ] * force_scaling),
                    ifelse(prescribed_forces[2], NaN, x[icol+1] * force_scaling),
                    ifelse(prescribed_forces[3], NaN, x[icol+2] * force_scaling))
        M_p[ip] = SVector(ifelse(prescribed_forces[4], NaN, x[icol+3] * force_scaling),
                    ifelse(prescribed_forces[5], NaN, x[icol+4] * force_scaling),
                    ifelse(prescribed_forces[6], NaN, x[icol+5] * force_scaling))

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(theta_p[ip])
        theta_p[ip] *= scaling

    end

    # element state variables
    u_e = Vector{SVector{3,TF}}(undef, ne)
    theta_e = Vector{SVector{3,TF}}(undef, ne)
    F_e = Vector{SVector{3,TF}}(undef, ne)
    M_e = Vector{SVector{3,TF}}(undef, ne)
    V_e = Vector{SVector{3,TF}}(undef, ne)
    Omega_e = Vector{SVector{3,TF}}(undef, ne)

    for ie = 1:ne

        icol = icol_elem[ie]

        u_e[ie] = SVector(x[icol], x[icol+1], x[icol+2])
        theta_e[ie] = SVector(x[icol+3], x[icol+4], x[icol+5])
        F_e[ie] = SVector(x[icol+6], x[icol+7], x[icol+8]) .* force_scaling
        M_e[ie] = SVector(x[icol+9], x[icol+10], x[icol+11]) .* force_scaling
        V_e[ie] = SVector(x[icol+12], x[icol+13], x[icol+14])
        Omega_e[ie] = SVector(x[icol+15], x[icol+16], x[icol+17])

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(theta_e[ie])
        theta_e[ie] *= scaling

    end

    return (u_p = u_p, theta_p = theta_p, F_p = F_p, M_p = M_p,
        u_e = u_e, theta_e = theta_e, F_e = F_e, M_e = M_e, V_e = V_e, Omega_e = Omega_e)
end

# convenience function for separating this model's input vector
function gxbeam_separate_inputs(y, displacement, icol_point, icol_elem)

    TF = eltype(y)

    np = length(icol_point)
    ne = length(icol_elem)

    # separate prescribed conditions
    u_p = Vector{SVector{3,TF}}(undef, np)
    theta_p = Vector{SVector{3,TF}}(undef, np)
    F_p = Vector{SVector{3,TF}}(undef, np)
    M_p = Vector{SVector{3,TF}}(undef, np)
    for ip = 1:np
        prescribed_forces = SVector{6}(view(displacement, :, ip)) .== false
        F_p[ip] = SVector(
            ifelse(prescribed_forces[1], y[6*(ip-1) + 1], NaN),
            ifelse(prescribed_forces[2], y[6*(ip-1) + 2], NaN),
            ifelse(prescribed_forces[3], y[6*(ip-1) + 3], NaN))
        M_p[ip] = SVector(
            ifelse(prescribed_forces[4], y[6*(ip-1) + 4], NaN),
            ifelse(prescribed_forces[5], y[6*(ip-1) + 5], NaN),
            ifelse(prescribed_forces[6], y[6*(ip-1) + 6], NaN))
        u_p[ip] = SVector(
            ifelse(prescribed_forces[1], NaN, y[6*(ip-1) + 1]),
            ifelse(prescribed_forces[2], NaN, y[6*(ip-1) + 2]),
            ifelse(prescribed_forces[3], NaN, y[6*(ip-1) + 3]))
        theta_p[ip] = SVector(
            ifelse(prescribed_forces[4], NaN, y[6*(ip-1) + 4]),
            ifelse(prescribed_forces[5], NaN, y[6*(ip-1) + 5]),
            ifelse(prescribed_forces[6], NaN, y[6*(ip-1) + 6]))
    end

    # separate distributed loads
    f_d = Vector{SVector{3,TF}}(undef, ne)
    m_d = Vector{SVector{3,TF}}(undef, ne)
    for ie = 1:ne
        f_d[ie] = SVector(
            y[6*np + 6*(ie-1) + 1], 
            y[6*np + 6*(ie-1) + 2], 
            y[6*np + 6*(ie-1) + 3])
        m_d[ie] = SVector(
            y[6*np + 6*(ie-1) + 4], 
            y[6*np + 6*(ie-1) + 5], 
            y[6*np + 6*(ie-1) + 6])
    end

    # separate point masses
    m_m = Vector{TF}(undef, ne)
    p_m = Vector{SVector{3,TF}}(undef, ne)
    I_m = Vector{SMatrix{3,3,TF}}(undef, ne)
    for ie = 1:ne
        m_m[ie] = y[6*np + 6*ne + 10*(ie-1) + 1]
        p_m[ie] = SVector(
            y[6*np + 6*ne + 10*(ie-1) + 2], 
            y[6*np + 6*ne + 10*(ie-1) + 3], 
            y[6*np + 6*ne + 10*(ie-1) + 4])
        I11 = y[6*np + 6*ne + 10*(ie-1) + 5]
        I22 = y[6*np + 6*ne + 10*(ie-1) + 6]
        I33 = y[6*np + 6*ne + 10*(ie-1) + 7]
        I12 = y[6*np + 6*ne + 10*(ie-1) + 8]
        I13 = y[6*np + 6*ne + 10*(ie-1) + 9]
        I23 = y[6*np + 6*ne + 10*(ie-1) + 10]
        I_m[ie] = @SMatrix [I11 I12 I13; I12 I22 I23; I13 I23 I33]
    end

    # gravity
    gvec =  view(y, 6*np + 6*ne + 10*ne +  1 : 6*np + 6*ne + 10*ne +  3)
    
    # linear/angular velocity
    v =     view(y, 6*np + 6*ne + 10*ne +  4 : 6*np + 6*ne + 10*ne +  6)
    omega = view(y, 6*np + 6*ne + 10*ne +  7 : 6*np + 6*ne + 10*ne +  9)
    
    # linear/angular acceleration
    a =     view(y, 6*np + 6*ne + 10*ne + 10 : 6*np + 6*ne + 10*ne + 12)
    alpha = view(y, 6*np + 6*ne + 10*ne + 13 : 6*np + 6*ne + 10*ne + 15)

    return (u_p=u_p, theta_p=theta_p, F_p=F_p, M_p=M_p, 
        f_d=f_d, m_d=m_d, m_m=m_m, p_m=p_m, I_m=I_m, 
        gvec=gvec, v=v, omega=omega, a=a, alpha=alpha)
end

# convenience function for separating this model's parameter vector
function gxbeam_separate_parameters(p, start, stop, icol_point, icol_elem)
    np = length(icol_point)
    ne = length(icol_elem)
    return (assembly = gxbeam_assembly(p, np, ne, start, stop),)
end

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