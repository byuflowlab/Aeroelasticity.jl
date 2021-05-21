"""
    GeometricallyExactBeam <: StructuralModel

Geometrically exact beam theory model, as implemented by the GXBeam package.
"""
struct GeometricallyExactBeam{TF, TV, TM} <: StructuralModel
    assembly::GXBeam.Assembly{TF}
    system::GXBeam.System{TF, TV, TM}
end
isinplace(::GeometricallyExactBeam) = true
islinear(::GeometricallyExactBeam) = false
is2D(::GeometricallyExactBeam) = false

# state rate terms form LHS
gebt_lhs!(Mdu, M, du) = mul!(Mdu, M, -du)

# steady state equations form RHS
function gebt_rhs!(Mdu, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    GXBeam.steady_state_system_residual!(Mdu, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return Mdu
end

# jacobian of RHS
function gebt_jacobian!(K, u, assembly, prescribed,
    distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    GXBeam.steady_state_system_jacobian!(K, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    return K
end

# mass matrix is calculated using GXBeam's internal function
function gebt_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
    irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    system_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    M .*= -1

    return M
end

function gebt_loads!(dq, r, system, assembly)

    # index in load array
    ir = 0

    for ielem = 1:nelem
        # get point indices
        irow_p1 = system.irow_pt[assembly.start[ielem]]
        irow_p2 = system.irow_pt[assembly.stop[ielem]]

        # extract loads on this element
        f = SVector(r[ir+1], r[ir+2], r[ir+3])
        m = SVector(r[ir+4], r[ir+5], r[ir+6])

        # get integrated distributed loads
        f1 = f2 = ΔL*f/2
        m1 = m2 = ΔL*m/2

        # subtract from residual vector
        dq[irow_p1:irow_p1+2] .-= f1 ./ system.force_scaling
        dq[irow_p1+3:irow_p1+5] .-= m1 ./ system.force_scaling
        dq[irow_p2:irow_p2+2] .-= f2 ./ system.force_scaling
        dq[irow_p2+3:irow_p2+5] .-= m2 ./ system.force_scaling

        # move to next set of distributed loads
        ir += 6
    end

    return dq
end

function gebt_load_jacobian(system, assembly)

    # index in load array
    ir = 0

    for ielem = 1:nelem
        # get point indices
        irow_p1 = system.irow_pt[assembly.start[ielem]]
        irow_p2 = system.irow_pt[assembly.stop[ielem]]

        # extract loads on this element
        f_r = Diagonal((@SVector ones(eltype(system), 3)))
        m_r = Diagonal((@SVector ones(eltype(system), 3)))

        # get integrated distributed loads
        f1_r = f2_r = ΔL*f_r/2
        m1_r = m2_r = ΔL*m_r/2

        # subtract from residual vector
        D[irow_p1:irow_p1+2, ir+1:ir+3] .= -f1_r ./ system.force_scaling
        D[irow_p1+3:irow_p1+5, ir+4:ir+6] .= -m1_r ./ system.force_scaling
        D[irow_p2:irow_p2+2, ir+1:ir+3] .= -f2_r ./ system.force_scaling
        D[irow_p2+3:irow_p2+5, ir+4:ir+6] .= -m2_r ./ system.force_scaling

        # move to next set of distributed loads
        ir += 6
    end

    return D
end

function coupled_rhs!(aero, stru::GeometricallyExactBeam, du, u, p, t)

    # extract GXBeam assembly and system
    assembly, system = stru.assembly, stru.system

    # number of aerodynamic and structural state variables
    Nq = length(system_state(system))
    Nλ = length(u) - Nq

    # extract states and state rates
    q, λ = view(u, 1:Nq), view(u, Nq+1:Nq+Nλ)
    dq, dλ = view(du, 1:Nq), view(du, Nq+1:Nq+Nλ)

    # extract structural parameters
    prescribed, distributed, x0, v0, ω0 = p

    # extract GXBeam control parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # get parameters for this point in time
    prescribed = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
    distributed = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
    x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
    v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
    ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

    # calculate mass matrix multiplied structural state rates
    gebt_rhs!(dq, q, assembly, prescribed, distributed, force_scaling,
        mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, x0, v0, ω0)

    # calculate aerodynamic loads
    get_loads!(aero, stru, r, u, p, t)

    # add contribution to structural state rates from aerodynamic loads
    gebt_loads!(dq, r, system, assembly)

    # calculate mass matrix multiplied aerodynamic state rates
    get_aerodynamic_rhs!(dλ, aero, stru, u, p, t)

    # return modified state rate vector
    return du
end

function coupled_jacobian!(aero, stru::GeometricallyExactBeam, J, u, p, t)

    # extract GXBeam assembly and system
    assembly, system = model.assembly, model.system

    # number of aerodynamic and structural state variables
    Nq = length(system_state(system))
    Nλ = length(u) - Nq

    # extract states and jacobians
    q, λ = view(u, 1:Nq), view(u, Nq+1:Nq+Nλ)
    Jss, Jsa = view(J, 1:Nq, 1:Nq), view(J, 1:Nq, Nq+1:Nq+Nλ)
    Jas, Jaa = view(J, Nq+1:Nq+Nλ, 1:Nq), view(J, Nq+1:Nq+Nλ, Nq+1:Nq+Nλ)

    # extract structural parameters
    prescribed, distributed, x0, v0, ω0 = p

    # extract GXBeam control parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # get parameters for this point in time
    prescribed = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
    distributed = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
    x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
    v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
    ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

    # calculate jacobian contributions due to the structure
    gebt_jacobian!(Jss, q, assembly, prescribed, distributed, force_scaling,
        mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam, x0, v0, ω0)

    # calculate jacobian contributions due to aerodynamic loads
    D = gebt_load_jacobian(system, assembly)
    Jls, Jla = get_load_jacobians(aero, stru, u, p, t)
    mul!(Jss, D, Jls, 1, 1)
    mul!(Jsa, D, Jla)

    # calculate jacobian components using the aerodynamic model
    get_aerodynamic_jacobians!(aero, stru, Jas, Jaa, u, p, t)

    # return modified jacobian
    return J
end

function coupled_mass_matrix!(aero, stru::GeometricallyExactBeam, M, u, p, t)
    # extract GXBeam assembly and system
    assembly, system = model.assembly, model.system

    # number of aerodynamic and structural state variables
    Nq = length(system_state(system))
    Nλ = length(u) - Nq

    # extract states and jacobians
    q, λ = view(u, 1:Nq), view(u, Nq+1:Nq+Nλ)
    Mss, Msa = view(M, 1:Nq, 1:Nq), view(M, 1:Nq, Nq+1:Nq+Nλ)
    Mas, Maa = view(M, Nq+1:Nq+Nλ, 1:Nq), view(M, Nq+1:Nq+Nλ, Nq+1:Nq+Nλ)

    # extract structural parameters
    prescribed, distributed, x0, v0, ω0 = p

    # extract GXBeam control parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # get parameters for this point in time
    prescribed = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
    distributed = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
    x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
    v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
    ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

    # add contributions to mass matrix components from the structural model
    gebt_mass_matrix!(Mss, q, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    # calculate mass matrix contributions due to aerodynamic loads
    D = gebt_load_jacobian(system, assembly)
    Mls, Mla = get_load_mass_matrices(aero, stru, u, p, t)
    mul!(Mss, D, Mls, 1, 1)
    mul!(Msa, D, Mla)

    # calculate mass matrix components using the aerodynamic model
    get_aerodynamic_jacobians!(aero, stru, Mas, Maa, u, p, t)

    # return the calculated mass matrix
    return M
end

# --- Structural Interface --- #

"""
    get_ode(model::GeometricallyExactBeam)

Return an in-place ODE function of the form `f!(du, u, p, t)` for a typical section
model with states `u = q` and parameters `p = [prescribed, distributed, x0, v0, ω0]`
"""
get_ode(model::GeometricallyExactBeam) = (du, u, p, t) -> get_structural_rates!(model, du, u, p, t)

function get_structural_rates!(model::GeometricallyExactBeam, du, u, p, t)
    # extract GXBeam assembly and system
    assembly, system = stru.assembly, stru.system
    # extract structural states
    q = view(u, 1:length(system_state(system)))
    # extract structural parameters
    prescribed, distributed, x0, v0, ω0 = p
    # extract GXBeam control parameters from the system
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    # get parameters for this point in time
    prescribed = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
    distributed = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
    x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
    v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
    ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
    # calculate mass matrix multiplied state rates
    gebt_rhs!(du, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
    return du
end

function get_mass_matrix!(model::GeometricallyExactBeam, M, u, p, t)
    gebt_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
    irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    # extract GXBeam assembly and system
    assembly, system = stru.assembly, stru.system
    # extract GXBeam control parameters from the system
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    gebt_mass_matrix!(M, u, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    return M
end

function get_jacobian!(model::GeometricallyExactBeam, J, u, p, t)
    # extract GXBeam assembly and system
    assembly, system = stru.assembly, stru.system
    # extract structural states
    q = view(u, 1:length(system_state(system)))
    # extract structural parameters
    prescribed, distributed, x0, v0, ω0 = p
    # extract GXBeam control parameters from the system
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    # get parameters for this point in time
    prescribed = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
    distributed = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
    x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
    v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
    ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
    # calculate jacobian of mass matrix multiplied state rates
    gebt_jacobian!(J, u, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
    return J
end

# TODO: Add get_parameter_jacobian

# TODO: Add get_syms
get_syms(::GeometricallyExactBeam) = nothing

# --- Coupled Interface --- #

"""
    get_ode(::TypicalSection)

Return an in-place ODE function of the form `f!(du, u, p, t)` for a typical section
model with states `u = q` and parameters `p = [prescribed, distributed, x0, v0, ω0]`
"""
function get_ode(aero::AerodynamicModel, stru::GeometricallyExactBeam)
    (du, u, p, t) -> get_aeroelastic_rates!(aero, stru, du, u, p, t)
end

function get_aeroelastic_rates!(aero::AerodynamicModel, stru::GeometricallyExactBeam, du, u, p, t)
    return coupled_rhs!(aero, stru, du, u, p, t)
end

function get_mass_matrix!(aero::AerodynamicModel, stru::GeometricallyExactBeam, M, u, p, t)
    return coupled_mass_matrix!(aero, stru, M, u, p, t)
end

function get_jacobian!(aero::AerodynamicModel, stru::GeometricallyExactBeam, J, u, p, t)
    return coupled_jacobian!(aero, stru, J, u, p, t)
end

# TODO: Add get_parameter_jacobian

get_syms(aero, ::GeometricallyExactBeam) = nothing
