"""
    GeometricallyExactBeam <: AbstractModel

Geometrically exact beam theory model, as implemented by the GXBeam package.
State variables are as defined by GXBeam. Inputs correspond to the (non-follower)
prescribed conditions on each point with prescribed conditions, followed by
constant (non-follower) distributed aerodynamic loads on each beam element with
distributed loads, followed by the model origin, linear velocity, and
angular velocity.  Parameters correspond to assembly point locations
``\\begin{bmatrix} x & y & z \\end{bmatrix}^T`` for each point in the assembly,
followed by the elements of the compliance matrix for each beam element in the
assembly ``\\begin{bmatrix} C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16} & C_{22}
& C_{23} & C_{24} & C_{25} & C_{26} & C_{33} & C_{34} & C_{35} & C_{36} & C_{44}
& C_{45} & C_{46} & C_{55} & C_{56} & C_{66} \end{bmatrix}^T``, followed by the
mass matrix parameters ``\\begin{bmatrix} μ & x_{m2} & x_{m3} & i_{22} & i_{33}
& i_{23} \\end{bmatrix}^T`` for each beam element in the assembly, followed by
the elements of the global to local frame transformation matrix
``\\begin{bmatrix} R_{11} & R{12} & R_{13} & R_{21} & R_{22} & R_{23} & R_{31}
& R_{32} & R_{33} \end{bmatrix}^T`` for each beam element in the assembly.
"""
struct GeometricallyExactBeam{TF, TV, TM} <: AbstractModel
    system::GXBeam.System{TF, TV, TM}
    assembly::GXBeam.Assembly{TF}
    prescribed_conditions::Dict{Int,GXBeam.PrescribedConditions{TF}}
    distributed_loads::Dict{Int,GXBeam.DistributedLoads{TF}}
end

# --- Constructors --- #

"""
    GeometricallyExactBeam(system, assembly, prescribed_conditions, distributed_loads)

Construct a GeometricallyExactBeam structural model.
"""
GeometricallyExactBeam(system, assembly, prescribed_conditions, distributed_loads)

# --- Traits --- #

number_of_states(model::GeometricallyExactBeam) = length(model.system.x)
function number_of_inputs(model::GeometricallyExactBeam)
    return 6*length(model.prescribed) + 6*length(model.distributed) + 9
end
function number_of_parameters(model::GeometricallyExactBeam)
    return 3*length(model.assembly.points) + 36*length(model.assembly.elements)
end
inplaceness(::Type{GeometricallyExactBeam}) = InPlace()
mass_matrix_type(::Type{GeometricallyExactBeam}) = Varying()
state_jacobian_type(::Type{GeometricallyExactBeam}) = Varying()
input_jacobian_type(::Type{GeometricallyExactBeam}) = Varying()
input_dependence_type(::Type{GeometricallyExactBeam}) = Linear()

# --- Methods --- #

function get_mass_matrix!(M, model::GeometricallyExactBeam, q, r, p, t)

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

    # update assembly using the parameter vector
    if eltype(p) <: eltype(assembly)
        # update existing assembly
        update_assembly!(assembly, p)
    else
        # create assembly with updated type
        assembly = update_assembly(assembly, p)
    end

    # return mass matrix
    return gebt_mass_matrix!(M, q, assembly, force_scaling, mass_scaling,
        irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
end

function get_rates!(dq, model::GeometricallyExactBeam, q, r, p, t)

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

    # update assembly using parameter vector
    if eltype(p) <: eltype(assembly)
        # update existing assembly using the parameters
        update_assembly!(assembly, p)
    else
        # create new assembly from template assembly and the parameters
        assembly = update_assembly(assembly, p)
    end

    # get number of (non-follower) prescribed inputs
    nrp = 6*length(prescribed)

    # extract portion of inputs for prescribed conditions
    rp = view(r, 1:nrp)

    # update prescribed conditions using input vector
    if eltype(eltype(prescribed)) <: eltype(rp)
        update_prescribed!(prescribed, rp)
    else
        prescribed = update_prescribed(prescribed, rp)
    end

    # get number of (non-follower) distributed inputs
    nrd = 6*length(prescribed)

    # extract portion of inputs for distributed loads
    rd = view(r, nrp+1:nrp+nrd)

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(rd)
        update_distributed!(distributed, rd)
    else
        distributed = update_distributed(distributed, rd)
    end

    # set origin, linear velocity, and angular velocity from the input vector
    x0 = SVector(r[end-8], r[end-7], r[end-6])
    v0 = SVector(r[end-5], r[end-4], r[end-3])
    ω0 = SVector(r[end-2], r[end-1], r[end])

    # return mass matrix multiplied state rates
    return gxbeam_rates!(Mdu, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_state_jacobian!(J, model::GeometricallyExactBeam, q, r, p, t)

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

    # update assembly using parameter vector
    if eltype(p) <: eltype(assembly)
        # update existing assembly using the parameters
        update_assembly!(assembly, p)
    else
        # create new assembly from template assembly and the parameters
        assembly = update_assembly(assembly, p)
    end

    # get number of (non-follower) prescribed inputs
    nrp = 6*length(prescribed)

    # extract portion of inputs for prescribed conditions
    rp = view(r, 1:nrp)

    # update prescribed conditions using input vector
    if eltype(eltype(prescribed)) <: eltype(rp)
        update_prescribed!(prescribed, rp)
    else
        prescribed = update_prescribed(prescribed, rp)
    end

    # get number of (non-follower) distributed inputs
    nrd = 6*length(prescribed)

    # extract portion of inputs for distributed loads
    rd = view(r, nrp+1:nrp+nrd)

    # update distributed loads using input vector
    if eltype(eltype(distributed)) <: eltype(rd)
        update_distributed!(distributed, rd)
    else
        distributed = update_distributed(distributed, rd)
    end

    # set origin, linear velocity, and angular velocity from the input vector
    x0 = SVector(r[end-8], r[end-7], r[end-6])
    v0 = SVector(r[end-5], r[end-4], r[end-3])
    ω0 = SVector(r[end-2], r[end-1], r[end])

    # return jacobian of right hand side with respect to state variables
    return gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
end

function get_input_jacobian!(Jr, model::GeometricallyExactBeam, q, r, p, t)

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
    return gxbeam_input_jacobian(q, r, start, stop, prescribed_conditions,
        distributed_loads, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, force_scaling, mass_scaling)
end

# TODO: Add parameter jacobian

# --- Internal --- #

# jacobian of LHS wrt state rates
function gebt_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
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
function gxbeam_input_jacobian(q, r, start, stop, prescribed_conditions,
    distributed_loads, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
    icol_beam, force_scaling, mass_scaling)

    f! = (y, x) -> gxbeam_jac_vec!(y, x, r, start, stop, prescribed_conditions,
        distributed_loads,irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, force_scaling, mass_scaling)
    M = length(q)
    N = length(r)

    return LinearMap(f!, M, N; ismutating = true)
end

# TODO: Parameter jacobian

# input jacobian vector product
function gxbeam_jac_vec!(y, x, r, start, stop, prescribed_conditions, distributed_loads,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    force_scaling, mass_scaling)

    # initialize output array
    y .= 0

    # current input array index
    ir = 0

    # get all points with prescribed conditions
    prescribed_points = keys(prescribed_conditions)

    # get all elements with distributed loads
    distributed_elements = keys(distributed_loads)

    # loop through each prescribed conditions (in order)
    for ipoint = 1:npoint
        # check if a prescribed condition exists for this point
        if ipoint in prescribed_points
            # check whether a force is prescribed
            prescribed_force = prescribed_conditions[ipoint].force
            # search for beams that are connected to this point
            for ibeam = 1:nbeam
                # check left side of beam
                if ipoint == beam_start[ibeam]
                    # add to residual equations for the beam endpoint
                    side = -1
                    irow_b = irow_beam1[ibeam]
                    if irow_b == irow_p
                        # add jacobian entries for equilibrium and compatability equations
                        for i = 1:3 # forces/displacements
                            if prescribed_force[i]
                                # F[i] is prescribed, u[i] is a state variable
                                y[irow_b+i-1] += -x[ir+i-1] # F_F = I
                            else
                                # u[i] is prescribed, F[i] is a state variable
                                y[irow_b+i+5] += side*x[ir+i-1] # u_u = I
                            end
                        end
                        for i = 4:6 # moments/rotations
                            if prescribed_force[i]
                                # M[i-3] is prescribed, θ[i-3] is a state variable
                                y[irow_b+i-1] += -x[ir+i-1] # M_M = I
                            else
                                # θ[i-3] is prescribed, M[i-3] is a state variable
                                y[irow_b+i+5] += side*x[ir+i-1] # θ_θ = I
                            end
                        end
                    else
                        # add jacobian entries for compatability equations only
                        for i = 1:3
                            if !prescribed_force[i]
                                # u[i] is prescribed, F[i] is a state variable
                                y[irow_b+i-1] += side*x[ir+i-1] # u_u = I
                            end
                        end
                        for i = 4:6
                            if !prescribed_force[i]
                                # θ[i-3] is prescribed, M[i-3] is a state variable
                                y[irow_b+i-1] += side*x[ir+i-1] # θ_θ = I
                            end
                        end
                    end
                end
                # check right side of beam
                if ipoint == beam_stop[ibeam]
                    # add to residual equations for the beam endpoint
                    side = 1
                    irow_b = irow_beam2[ibeam]
                    if irow_b == irow_p
                        # add jacobian entries for equilibrium and compatability equations
                        for i = 1:3 # forces/displacements
                            if prescribed_force[i]
                                # F[i] is prescribed, u[i] is a state variable
                                y[irow_b+i-1] += -x[ir+i-1] # F_F = I
                            else
                                # u[i] is prescribed, F[i] is a state variable
                                y[irow_b+i+5] += side*x[ir+i-1] # u_u = I
                            end
                        end
                        for i = 4:6 # moments/rotations
                            if prescribed_force[i]
                                # M[i-3] is prescribed, θ[i-3] is a state variable
                                y[irow_b+i-1] += -x[ir+i-1] # M_M = I
                            else
                                # θ[i-3] is prescribed, M[i-3] is a state variable
                                y[irow_b+i+5] += side*x[ir+i-1] # θ_θ = I
                            end
                        end
                    else
                        # add jacobian entries for compatability equations only
                        for i = 1:3
                            if !prescribed_force[i]
                                # u[i] is prescribed, F[i] is a state variable
                                y[irow_b+i-1] += side*x[ir+i-1] # u_u = I
                            end
                        end
                        for i = 4:6
                            if !prescribed_force[i]
                                # θ[i-3] is prescribed, M[i-3] is a state variable
                                y[irow_b+i-1] += side*x[ir+i-1] # θ_θ = I
                            end
                        end
                    end
                end
            end
            # move to next set of prescribed conditions
            ir += 6
        end
    end

    # distributed loads
    for ielem = 1:nelem
        # get point indices
        irow_p1 = irow_pt[beam_start[ielem]]
        irow_p2 = irow_pt[beam_stop[ielem]]
        # initialize element load jacobians
        f_r = Diagonal((@SVector ones(3)))
        m_r = Diagonal((@SVector ones(3)))
        # jacobians after integration
        f1_r = f2_r = ΔL*f_r/2
        m1_r = m2_r = ΔL*m_r/2
        # insert jacobian contributions
        jacob[irow_p1:irow_p1+2, ir+1:ir+3] .= -f1_r ./ force_scaling
        jacob[irow_p1+3:irow_p1+5, ir+4:ir+6] .= -m1_r ./ force_scaling
        jacob[irow_p2:irow_p2+2, ir+1:ir+3] .= -f2_r ./ force_scaling
        jacob[irow_p2+3:irow_p2+5, ir+4:ir+6] .= -m2_r ./ force_scaling
        # move to next set of distributed loads
        ir += 6
    end

    return y
end

# update assembly using parameters
function update_assembly(assembly, p)
    new_points = SVector{3,eltype(p)}.(assembly.points)
    new_start = assembly.start
    new_stop = assembly.stop
    new_elements = Element{TF}.(assembly.elements)
    new_assembly = Assembly(new_points, new_start, new_stop, new_elements)
    update_assembly!(new_assembly, p)
    return new_assembly
end

# update assembly using parameters
function update_assembly!(assembly, p)
    # unpack assembly
    points = assembly.points
    elements = assembly.elements

    # number of points and elements
    npoint = length(points)
    nelem = length(elements)

    # number of entries in parameter vector corresponding to points
    npp = 3*npoint

    # number of entries in parameter vector corresponding to compliance matrices
    npc = 21*nelem

    # number of entries in parameter vector corresponding to inertial parameters
    npm = 6*nelem

    # number of entries in parameter vector corresponding to rotation matrices
    npr = 9*nelem

    # update each point
    li = LinearIndices(3, length(points))
    for (ipoint, point) in enumerate(points)
        # point indices in `r`
        ix = li[1,ipoint]
        iy = li[2,ipoint]
        iz = li[3,ipoint]
        # create new point
        new_point = p[SVector(ix, iy, iz)]
        # replace old point coordinates
        if ismutable(point)
            point .= new_point
        else
            points[ipoint] = new_point
        end
    end

    # update each element
    for (ielem, element) in enumerate(elements)
        # calculate length from point locations
        ΔL = norm(points[stop[ielem]] - points[start[ielem]])
        # calculate midpoint from point locations
        x = (points[stop[ielem]] + points[start[ielem]])/2
        # extract stiffness parameters
        offset = npp + 21*(ielem - 1)
        C_11 = p[offset+1]
        C_12 = p[offset+2]
        C_13 = p[offset+3]
        C_14 = p[offset+4]
        C_15 = p[offset+5]
        C_16 = p[offset+6]
        C_22 = p[offset+7]
        C_23 = p[offset+8]
        C_24 = p[offset+9]
        C_25 = p[offset+10]
        C_26 = p[offset+11]
        C_33 = p[offset+12]
        C_34 = p[offset+13]
        C_35 = p[offset+14]
        C_36 = p[offset+15]
        C_44 = p[offset+16]
        C_45 = p[offset+17]
        C_46 = p[offset+18]
        C_55 = p[offset+19]
        C_56 = p[offset+20]
        C_66 = p[offset+21]
        # construct compliance matrix
        C = @SMatrix [C_11 C_12 C_13 C_14 C_15 C_16;
            C_12 C_22 C_23 C_24 C_25 C_26;
            C_13 C_23 C_33 C_34 C_35 C_36;
            C_14 C_24 C_34 C_44 C_45 C_46;
            C_15 C_25 C_35 C_45 C_55 C_56;
            C_16 C_26 C_36 C_46 C_56 C_66]
        # extract mass parameters
        offset = npp + npc + 6*(ielem - 1)
        μ & x_{m2} & x_{m3} & i_{22} & i_{33} & i_{23}
        μ = p[offset+1]
        x_m2 = p[offset+2]
        x_m3 = p[offset+3]
        i_22 = p[offset+4]
        i_33 = p[offset+5]
        i_23 = p[offset+6]
        # construct inverse mass matrix
        mass = @SMatrix [     μ       0      0           0 μ*x_m3 -μ*x_m2;
                              0       μ      0     -μ*x_m3      0       0;
                              0       0      μ      μ*x_m2      0       0;
                              0 -μ*x_m3 μ*x_m2 i_22 + i_33      0       0;
                         μ*x_m3       0      0           0   i_22   -i_23;
                        -μ*x_m2       0      0           0  -i_23    i_33]
        minv = inv(mass)
        # extract transformation matrix parameters
        offset = npp + npc + npm + 9*(ielem - 1)
        R_11 = p[offset+1]
        R_21 = p[offset+2]
        R_31 = p[offset+3]
        R_12 = p[offset+4]
        R_22 = p[offset+5]
        R_32 = p[offset+6]
        R_31 = p[offset+7]
        R_32 = p[offset+8]
        R_33 = p[offset+9]
        # construct transformation matrix
        Cab = @SMatrix [R_11 R_12 R_13; R_21 R_22 R_23; R_31 R_32 R_33]
        # replace old element with new element
        elements[ielem] = Element(L, x, C, minv, Cab)
    end

    return assembly
end

# update prescribed conditions using inputs
function update_prescribed(prescribed, r::AbstractVector{T}) where T
    new_prescribed = Dict{Int, PrescribedConditions{T}}()
    for (key, value) in prescribed
        new_prescribed[key] = value
    end
    update_prescribed!(prescribed, r)
    return new_prescribed
end

# update prescribed conditions using inputs
function update_prescribed!(prescribed, r)
    prescribed_points = keys(prescribed)
    largest_point = maximum(prescribed_points)
    ir = 0
    for ipoint in 1:largest_point
        if ipoint in prescribed_points
            force = prescribed.force
            value = view(r, ir+1:ir+6)
            follower = prescribed.follower
            prescribed[ipoint] = PrescribedConditions(force, value, follower)
            ir += 6
        end
    end
    return prescribed
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

# function gebt_loads!(dq, r, system, assembly)
#
#     # index in load array
#     ir = 0
#
#     for ielem = 1:nelem
#         # get point indices
#         irow_p1 = system.irow_pt[assembly.start[ielem]]
#         irow_p2 = system.irow_pt[assembly.stop[ielem]]
#
#         # extract loads on this element
#         f = SVector(r[ir+1], r[ir+2], r[ir+3])
#         m = SVector(r[ir+4], r[ir+5], r[ir+6])
#
#         # get integrated distributed loads
#         f1 = f2 = ΔL*f/2
#         m1 = m2 = ΔL*m/2
#
#         # subtract from residual vector
#         dq[irow_p1:irow_p1+2] .-= f1 ./ system.force_scaling
#         dq[irow_p1+3:irow_p1+5] .-= m1 ./ system.force_scaling
#         dq[irow_p2:irow_p2+2] .-= f2 ./ system.force_scaling
#         dq[irow_p2+3:irow_p2+5] .-= m2 ./ system.force_scaling
#
#         # move to next set of distributed loads
#         ir += 6
#     end
#
#     return dq
# end
