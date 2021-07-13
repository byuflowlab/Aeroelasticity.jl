"""
    LiftingLine{N,T} <: AbstractModel

Lifting line model with `N` cross sections, using the aerodynamic models in `T`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated
"""
struct LiftingLine{N,T} <: AbstractModel
    models::T
end

"""
    LiftingLineSection <: AbstractModel

Lifting line section model with state variables ``q = \\begin{bmatrix} v_x & v_y &
v_z & \\omega_x & \\omega_y & \\omega_z \\end{bmatrix}^T``, inputs
``r = \\begin{bmatrix} L' & Y' & D' & Mx, My, Mz \\end{bmatrix}^T``, and no
parameters.  Two-dimensional aerodynamic models may be extended to
three-dimensional models by coupling with this model.  Note that
this model has no rate equations of its own since its state variables are
defined as functions of the 3D structural model's state variables.
"""
struct LiftingLineSection <: AbstractModel end

number_of_states(::Type{<:LiftingLineSection}) = 6
number_of_inputs(::Type{<:LiftingLineSection}) = 6
number_of_parameters(::Type{<:LiftingLineSection}) = 0
inplaceness(::Type{LiftingLineSection}) = OutOfPlace()

# --- Constructors --- #

"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of aerodynamic models.
"""
LiftingLine(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

LiftingLine{N}(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

"""
    LiftingLine{N}(model)

Construct a lifting line aerodynamic model using `N` instances of `model`.
"""
function LiftingLine{N}(model) where {N,T}
    models = ntuple(i->model, N)
    return LiftingLine(models)
end

# --- Traits --- #

number_of_states(model::LiftingLine) = sum(number_of_states.(model.models))
number_of_inputs(model::LiftingLine) = sum(number_of_inputs.(model.models))
number_of_parameters(model::LiftingLine) = sum(number_of_parameters.(model.models))
inplaceness(::Type{LiftingLine{N,T}}) where {N,T} = InPlace()

function mass_matrix_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(mass_matrix_type.(model_types)))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types)))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types)))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function input_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(input_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(input_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(input_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(input_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(input_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- Methods --- #

function get_rates!(du, model::LiftingLine, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dus = view.(Ref(du), iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_rates!.(dus, models, us, ys, ps, t)

    return du
end

function get_mass_matrix!(M, model::LiftingLine)

    M .= 0

    models = model.models

    iu = state_indices(model.models)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models)

    return M
end

function get_mass_matrix!(M, model::LiftingLine, u, y, p, t)

    M .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Ms = view.(Ref(M), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_mass_matrix!.(Ms, models, us, ys, ps, t)

    return M
end

# --- Performance Overloads --- #

function get_state_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    iu = state_indices(model.models)

    Js = view.(Ref(J), iu, iu)

    get_state_jacobian!.(Js, models)

    return J
end

function get_state_jacobian!(J, model::LiftingLine, u, y, p, t)

    J .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_state_jacobian!.(Js, models, us, ys, ps, t)

    return J
end

function get_input_jacobian(model::LiftingLine)

    f! = (y, x) -> input_jacobian_product!(y, x, model)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

function get_input_jacobian(model::LiftingLine, λ, d, p, t)

    f! = (y, x) -> input_jacobian_product!(y, x, model, λ, d, p, t)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

# --- Unit Testing Methods --- #

function get_lhs(model::LiftingLine, du, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dus = view.(Ref(du), iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    return vcat(get_lhs.(models, dus, us, ys, ps, t)...)
end

# --- Rigid Body Coupling --- #

"""
    couple_models(aero::LiftingLine, stru::RigidBody)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a rigid body model.  This model introduces the additional parameters
``\\begin{bmatrix} \\Delta L & p_b & b_x & b_y & b_z \\end{bmatrix}^T`` for each
lifting line element, followed by the freestream velocity components
``\\begin{bmatrix} V_x & V_y & V_z \\end{bmatrix}^T``, air density ``\\rho``,
gravitational constant ``g``, inertial properties ``m, Ixx, Iyy, Izz, Ixz, Ixy,
Iyz``, and constant applied forces/moments (in the body frame) ``Fx, Fy, Fz, Mx,
My, Mz``.

** When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction **
"""
couple_models(aero::LiftingLine, stru::RigidBody)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(::Type{<:LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    return 13*N+18
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::RigidBody, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = q
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # extract coupled system parameters
    Vinf = SVector(p[npa+nps+13*N+1], p[npa+nps+13*N+2], p[npa+nps+13*N+3])
    ρ = p[npa+nps+13*N+4]
    g = p[npa+nps+13*N+5]
    m = p[npa+nps+13*N+6]
    Ixx = p[npa+nps+13*N+7]
    Iyy = p[npa+nps+13*N+8]
    Izz = p[npa+nps+13*N+9]
    Ixz = p[npa+nps+13*N+10]
    Ixy = p[npa+nps+13*N+11]
    Iyz = p[npa+nps+13*N+12]
    Fx = p[npa+nps+13*N+13]
    Fy = p[npa+nps+13*N+14]
    Fz = p[npa+nps+13*N+15]
    Mx = p[npa+nps+13*N+16]
    My = p[npa+nps+13*N+17]
    Mz = p[npa+nps+13*N+18]

    # add freestream velocity due to translation
    Vinf -= V

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # loop through each lifting line section
    for i = 1:N

        # extract length, position, and orientation of this section
        ΔL = p[npa+nps+13*(i-1)+1]
        pb = SVector(p[npa+nps+13*(i-1)+2], p[npa+nps+13*(i-1)+3], p[npa+nps+13*(i-1)+4])
        bx = SVector(p[npa+nps+13*(i-1)+5], p[npa+nps+13*(i-1)+6], p[npa+nps+13*(i-1)+7])
        by = SVector(p[npa+nps+13*(i-1)+8], p[npa+nps+13*(i-1)+9], p[npa+nps+13*(i-1)+10])
        bz = SVector(p[npa+nps+13*(i-1)+11], p[npa+nps+13*(i-1)+12], p[npa+nps+13*(i-1)+13])

        # extract element aerodynamic states
        λi = SVector{Nλi[i]}(λs[i])

        # transformation matrix from the body frame to the local aerodynamic frame
        Cba = [bx by bz]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = Cba*Vinf

        # add velocity due to rotation
        vi -= Cba*cross(Ω, pb)

        # section angular velocity is body angular velocity
        ωi = Cba*Ω

        # calculate lifting line section inputs
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs(aero.models[i], LiftingLineSection(), ui, pi, t)

        # extract and save aerodynamic model inputs (in local aerodynamic frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments (in body frame)
        Fi = Cba'*ΔL*SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        Mi = Cba'*ΔL*SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # add contribution to total forces and moments
        Ftot += Fi
        Mtot += cross(pb, Fi) + Mi

    end

    # add graviational forces
    sϕ, cϕ = sincos(ϕr)
    sθ, cθ = sincos(θr)

    Fg = m*g*[-sθ, cθ*sϕ, cθ*cϕ]
    Mg = @SVector zeros(3)  # no gravitational moment

    # save rigid body inputs
    ys[1] = m
    ys[2] = Ixx
    ys[3] = Iyy
    ys[4] = Izz
    ys[5] = Ixz
    ys[6] = Ixy
    ys[7] = Iyz
    ys[8] = Ftot[1] + Fg[1] + Fx
    ys[9] = Ftot[2] + Fg[2] + Fy
    ys[10] = Ftot[3] + Fg[3] + Fz
    ys[11] = Mtot[1] + Mg[1] + Mx
    ys[12] = Mtot[2] + Mg[2] + My
    ys[13] = Mtot[3] + Mg[3] + Mz

    return y
end

function get_input_mass_matrix!(M, aero::LiftingLine{N,T}, stru::RigidBody,
    u, p, t) where {N,T}

    # start with zeroed out mass matrix
    M .= 0

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables and parameters
    λ = view(u, iλ) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = q
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # extract coupled system parameters
    Vinf = SVector(p[npa+nps+13*N+1], p[npa+nps+13*N+2], p[npa+nps+13*N+3])
    ρ = p[npa+nps+13*N+4]
    g = p[npa+nps+13*N+5]
    m = p[npa+nps+13*N+6]
    Ixx = p[npa+nps+13*N+7]
    Iyy = p[npa+nps+13*N+8]
    Izz = p[npa+nps+13*N+9]
    Ixz = p[npa+nps+13*N+10]
    Ixy = p[npa+nps+13*N+11]
    Iyz = p[npa+nps+13*N+12]
    Fx = p[npa+nps+13*N+13]
    Fy = p[npa+nps+13*N+14]
    Fz = p[npa+nps+13*N+15]
    Mx = p[npa+nps+13*N+16]
    My = p[npa+nps+13*N+17]
    Mz = p[npa+nps+13*N+18]

    # add freestream velocity due to translation
    Vinf -= V
    dVinf_dV = -SMatrix{3,3}(I)

    # initialize total forces and moments
    Ftot_dV = @SMatrix zeros(3, 3)
    Ftot_dΩ = @SMatrix zeros(3, 3)
    Mtot_dV = @SMatrix zeros(3, 3)
    Mtot_dΩ = @SMatrix zeros(3, 3)

    # loop through each lifting line / beam element
    for i = 1:N

        # extract length, position, and orientation of this section
        ΔL = p[npa+nps+13*(i-1)+1]
        pb = SVector(p[npa+nps+13*(i-1)+2], p[npa+nps+13*(i-1)+3], p[npa+nps+13*(i-1)+4])
        bx = SVector(p[npa+nps+13*(i-1)+5], p[npa+nps+13*(i-1)+6], p[npa+nps+13*(i-1)+7])
        by = SVector(p[npa+nps+13*(i-1)+8], p[npa+nps+13*(i-1)+9], p[npa+nps+13*(i-1)+10])
        bz = SVector(p[npa+nps+13*(i-1)+11], p[npa+nps+13*(i-1)+12], p[npa+nps+13*(i-1)+13])

        # extract element aerodynamic states
        λi = SVector{Nλi[i]}(λs[i])

        # transformation matrix from the body frame to the local aerodynamic frame
        Cba = [bx by bz]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = Cba*Vinf
        dv_dV = Cba*dVinf_dV

        # add velocity due to rotation
        vi -= Cba*cross(Ω, pb)
        dv_dΩ = Cba*GXBeam.tilde(pb)

        # section angular velocity is body angular velocity
        ωi = Cba*Ω
        dω_dΩ = Cba

        # calculate lifting line section mass matrices
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        Mi = get_input_mass_matrix(aero.models[i], LiftingLineSection(), ui, pi, t)

        # separate into component mass matrices
        d_dλ = SMatrix{Nyai[i], Nλi[i]}(view(Mi, 1:Nyai[i], 1:Nλi[i]))

        d_dv = SMatrix{Nyai[i], 3}(view(Mi, 1:Nyai[i], Nλi[i]+1 : Nλi[i]+3))
        d_dω = SMatrix{Nyai[i], 3}(view(Mi, 1:Nyai[i], Nλi[i]+4 : Nλi[i]+6))
        d_dV = d_dv * dv_dV
        d_dΩ = d_dv * dv_dΩ + d_dω * dω_dΩ

        f_dλ = SMatrix{3, Nλi[i]}(view(Mi, Nyai[i]+1 : Nyai[i]+3, 1:Nλi[i]))
        m_dλ = SMatrix{3, Nλi[i]}(view(Mi, Nyai[i]+4 : Nyai[i]+6, 1:Nλi[i]))

        f_dv = SMatrix{3, 3}(view(Mi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+1 : Nλi[i]+3))
        f_dω = SMatrix{3, 3}(view(Mi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+4 : Nλi[i]+6))
        f_dV = f_dv * dv_dV
        f_dΩ = f_dv * dv_dΩ + f_dω * dω_dΩ

        m_dv = SMatrix{3, 3}(view(Mi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+1 : Nλi[i]+3))
        m_dω = SMatrix{3, 3}(view(Mi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+4 : Nλi[i]+6))
        m_dV = m_dv * dv_dV
        m_dΩ = m_dv * dv_dΩ + m_dω * dω_dΩ

        # save aerodynamic input mass matrix entries (in local aerodynamic frame)
        M[iyas[i], iλs[i]] = d_dλ
        M[iyas[i], nλ+7:nλ+9] = d_dV
        M[iyas[i], nλ+10:nλ+12] = d_dΩ

        # extract local forces and moments (in body frame)
        Fi_dλ = Cba'*ΔL*f_dλ
        Fi_dV = Cba'*ΔL*f_dV
        Fi_dΩ = Cba'*ΔL*f_dΩ
        Mi_dλ = Cba'*ΔL*m_dλ
        Mi_dV = Cba'*ΔL*m_dV
        Mi_dΩ = Cba'*ΔL*m_dΩ

        # add contribution to total forces and moments
        M[nya+8 : nya+10, iλs[i]] = Fi_dλ
        Ftot_dV += Fi_dV
        Ftot_dΩ += Fi_dΩ

        M[nya+11 : nya+13, iλs[i]] = GXBeam.tilde(pb)*Fi_dλ + Mi_dλ
        Mtot_dV += GXBeam.tilde(pb)*Fi_dV + Mi_dV
        Mtot_dΩ += GXBeam.tilde(pb)*Fi_dΩ + Mi_dΩ

    end

    # save rigid body inputs
    M[nya+8 : nya+10, nλ+7:nλ+9] = Ftot_dV
    M[nya+8 : nya+10, nλ+10:nλ+12] = Ftot_dΩ
    M[nya+11 : nya+13, nλ+7:nλ+9] = Mtot_dV
    M[nya+11 : nya+13, nλ+10:nλ+12] = Mtot_dΩ

    return M
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::RigidBody,
    du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Nu = number_of_inputs(models)
    y = zeros(TF, Nu)

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dλ = view(du, iλ) # aerodynamic states rates
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dq = view(du, iq) # structural state rates
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    dλs = view.(Ref(dλ), iλs) # aerodynamic state rates for each section
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states and rates
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = q
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dq
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)
    dV = SVector(dur, dvr, dwr)
    dΩ = SVector(dpr, dqr, drr)

    # extract coupled system parameters
    Vinf = SVector(p[npa+nps+13*N+1], p[npa+nps+13*N+2], p[npa+nps+13*N+3])
    ρ = p[npa+nps+13*N+4]
    g = p[npa+nps+13*N+5]
    m = p[npa+nps+13*N+6]
    Ixx = p[npa+nps+13*N+7]
    Iyy = p[npa+nps+13*N+8]
    Izz = p[npa+nps+13*N+9]
    Ixz = p[npa+nps+13*N+10]
    Ixy = p[npa+nps+13*N+11]
    Iyz = p[npa+nps+13*N+12]
    Fx = p[npa+nps+13*N+13]
    Fy = p[npa+nps+13*N+14]
    Fz = p[npa+nps+13*N+15]
    Mx = p[npa+nps+13*N+16]
    My = p[npa+nps+13*N+17]
    Mz = p[npa+nps+13*N+18]

    # add freestream velocity due to translation
    Vinf -= V

    # set freestream acceleration
    dVinf = -dV

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # loop through each lifting line / beam element
    for i = 1:N

        # extract length, position, and orientation of this section
        ΔL = p[npa+nps+13*(i-1)+1]
        pb = SVector(p[npa+nps+13*(i-1)+2], p[npa+nps+13*(i-1)+3], p[npa+nps+13*(i-1)+4])
        bx = SVector(p[npa+nps+13*(i-1)+5], p[npa+nps+13*(i-1)+6], p[npa+nps+13*(i-1)+7])
        by = SVector(p[npa+nps+13*(i-1)+8], p[npa+nps+13*(i-1)+9], p[npa+nps+13*(i-1)+10])
        bz = SVector(p[npa+nps+13*(i-1)+11], p[npa+nps+13*(i-1)+12], p[npa+nps+13*(i-1)+13])

        # extract element aerodynamic states and rates
        λi = SVector{Nλi[i]}(λs[i])
        dλi = SVector{Nλi[i]}(dλs[i])

        # transformation matrix from the body frame to the local aerodynamic frame
        Cba = [bx by bz]

        # transform freestream velocity and acceleration from body frame to local aerodynamic frame
        vi = Cba * Vinf
        dvi = Cba * dVinf

        # add velocity due to rotation
        vi -= Cba*cross(Ω, pb)
        dvi -= Cba*cross(dΩ, pb)

        # section angular velocity is body angular velocity
        ωi = Cba*Ω
        dωi = Cba*dΩ

        # calculate lifting line section inputs
        dui = vcat(dλi, dvi, dωi) # section state rates
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs_from_state_rates(aero.models[i], LiftingLineSection(),
            dui, ui, pi, t)

        # extract and save aerodynamic model inputs (in local aerodynamic frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments (in body frame)
        Fi = Cba'*ΔL*SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        Mi = Cba'*ΔL*SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # add contribution to total forces and moments
        Ftot += Fi
        Mtot += cross(pb, Fi) + Mi

    end

    # save total forces and moments
    ys[8:10] = Ftot
    ys[11:13] = Mtot

    return y
end

# --- Geometrically Exact Beam Theory Coupling --- #

"""
    couple_models(aero::LiftingLine, stru::GEBT)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory model.  This model introduces additional
parameters corresponding to the freestream velocity components ``\\begin{bmatrix}
V_x & V_y & V_z \\end{bmatrix}^T``, freestream air density ``\\rho``, external
forces ``F_{x,i}``, ``F_{y,i}``, ``F_{z,i}``, ``M_{x,i}``, ``M_{y,i}``,
``M_{z,i}`` or displacements ``u_{x,i}``, ``u_{y,i}``, ``u_{z,i}``,
``\\theta_{x,i}``, ``\\theta_{y,i}``, ``\\theta_{z,i}`` applied to each node,
constant distributed loads ``f_{x,i}``, ``f_{y,i}``, ``f_{z,i}``, ``m_{x,i}``,
``m_{y,i}``, ``m_{z,i}`` applied to each beam element, and the linear and
angular velocity ``V_x``, ``V_y``, ``V_z``, ``\\Omega_x``, ``\\Omega_y``,
``\\Omega_z`` of the system.

** When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface **
"""
couple_models(aero::LiftingLine, stru::GEBT)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(aero::LiftingLine, stru::GEBT)
    return 4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 6
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])
    ρ = p[npa+nps+4]

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from structural parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start,
        stru.element_stop)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = 6*(ip-1)
        poff = npa + nps + 4 + 6*(ip-1)
        ys[yoff+1] = p[poff+1]
        ys[yoff+2] = p[poff+2]
        ys[yoff+3] = p[poff+3]
        ys[yoff+4] = p[poff+4]
        ys[yoff+5] = p[poff+5]
        ys[yoff+6] = p[poff+6]
    end

    # save lifting line/beam element loads/displacements
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic states
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add local freestream linear velocities due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # calculate lifting line section inputs
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs(aero.models[i], LiftingLineSection(), ui, pi, t)

        # save aerodynamic inputs for this element
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # constant distributed loads (in body frame)
        poff = npa + nps + 4 + 6*npoint + 6*(i-1)
        fic = SVector(p[poff+1], p[poff+2], p[poff+3])
        mic = SVector(p[poff+4], p[poff+5], p[poff+6])

        # aerodynamic distributed loads (in body frame)
        fia = CtCab*R'*SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mia = CtCab*R'*SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        ys[yoff+1 : yoff+3] = fic + fia
        ys[yoff+4 : yoff+6] = mic + mia

    end

    # save body frame linear/angular velocities
    yoff = 6*npoint + 6*nelem
    poff = npa + nps + 4 + 6*npoint + 6*nelem
    ys[yoff+1] = p[poff+1]
    ys[yoff+2] = p[poff+2]
    ys[yoff+3] = p[poff+3]
    ys[yoff+4] = p[poff+4]
    ys[yoff+5] = p[poff+5]
    ys[yoff+6] = p[poff+6]

    return y
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # start with zero valued mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])
    ρ = p[npa+nps+4]

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start, stru.element_stop)

    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract aerodynamic state variables corresponding to this element
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add local freestream linear velocities due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        dv_dP = -R * element.minv11 * stru.mass_scaling
        dv_dH = -R * element.minv12 * stru.mass_scaling

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        dω_dP = R * element.minv12' * stru.mass_scaling
        dω_dH = R * element.minv22 * stru.mass_scaling

        # calculate lifting line section input mass matrix
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        Myi = get_input_mass_matrix(aero.models[i], LiftingLineSection(), ui, pi, t)

        # separate into component mass matrices
        d_dλ = SMatrix{Nyai[i], Nλi[i]}(view(Myi, 1:Nyai[i], 1:Nλi[i]))

        d_dv = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+1 : Nλi[i]+3))
        d_dω = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+4 : Nλi[i]+6))
        d_dP = d_dv * dv_dP + d_dω * dω_dP
        d_dH = d_dv * dv_dH + d_dω * dω_dH

        f_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+1 : Nyai[i]+3, 1:Nλi[i]))
        m_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+4 : Nyai[i]+6, 1:Nλi[i]))

        f_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+1 : Nλi[i]+3))
        f_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+4 : Nλi[i]+6))
        f_dP = f_dv * dv_dP + f_dω * dω_dP
        f_dH = f_dv * dv_dH + f_dω * dω_dH

        m_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+1 : Nλi[i]+3))
        m_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+4 : Nλi[i]+6))
        m_dP = m_dv * dv_dP + m_dω * dω_dP
        m_dH = m_dv * dv_dH + m_dω * dω_dH

        # save aerodynamic input mass matrix entries (in local aerodynamic frame)
        icol = stru.icol_beam[i]
        My[iyas[i], iλs[i]] = d_dλ
        My[iyas[i], nλ+icol+12:nλ+icol+14] = d_dP
        My[iyas[i], nλ+icol+15:nλ+icol+17] = d_dH

        # save load mass matrix entries (in body frame)
        offset = nya + 6*(i-1)
        My[offset+1 : offset+3, iλs[i]] = CtCab*R'*f_dλ
        My[offset+4 : offset+6, iλs[i]] = CtCab*R'*m_dλ
        My[offset+1 : offset+3, nλ+icol+12:nλ+icol+14] = CtCab*R'*f_dP
        My[offset+1 : offset+3, nλ+icol+15:nλ+icol+17] = CtCab*R'*f_dH
        My[offset+4 : offset+6, nλ+icol+12:nλ+icol+14] = CtCab*R'*m_dP
        My[offset+4 : offset+6, nλ+icol+15:nλ+icol+17] = CtCab*R'*m_dH

    end

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::GEBT, du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Nu = number_of_inputs(models)
    y = zeros(TF, Nu)

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dλ = view(du, iλ) # aerodynamic states rates
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dq = view(du, iq) # structural state rates
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    dλs = view.(Ref(dλ), iλs) # aerodynamic state rates for each section
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])
    ρ = p[npa+nps+4]

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start, stru.element_stop)


    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic state rates and states
        dλi = SVector{Nλi[i]}(dλs[i])
        λi = SVector{Nλi[i]}(λs[i])

        # get element structural state rates and states
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        du_elem = SVector(dq[icol], dq[icol+1], dq[icol+2])
        dθ_elem = SVector(dq[icol+3], dq[icol+4], dq[icol+5])
        dF_elem = SVector(dq[icol+6], dq[icol+7], dq[icol+8]) .* stru.force_scaling
        dM_elem = SVector(dq[icol+9], dq[icol+10], dq[icol+11]) .* stru.force_scaling
        dP_elem = SVector(dq[icol+12], dq[icol+13], dq[icol+14]) .* stru.mass_scaling
        dH_elem = SVector(dq[icol+15], dq[icol+16], dq[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling
        dθ_elem *= scaling

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add velocity due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        # calculate freestream acceleration due to surface motion
        dvi = -R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # calculate local section angular accelerations due to surface motion
        dωi = R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # calculate lifting line section inputs from state rate contributions
        dui = vcat(dλi, dvi, dωi) # section state rates
        ui = vcat(λi, vi, ωi) # section states
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs_from_state_rates(aero.models[i], LiftingLineSection(),
            dui, ui, pi, t)

        # extract and save aerodynamic model inputs (in local aerodynamic frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments per unit length (in local aerodynamic frame)
        fi = SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mi = SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # save force and moment per unit length (in body frame)
        ys[6*(i-1)+1 : 6*(i-1)+3] = CtCab*R'*fi
        ys[6*(i-1)+4 : 6*(i-1)+6] = CtCab*R'*mi

    end

    return y
end

"""
    couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory and rigid body model.  This model introduces additional
parameters corresponding to the freestream velocity components ``\\begin{bmatrix}
V_x & V_y & V_z \\end{bmatrix}^T``, air density ``\\rho``, gravitational
constant ``g``, and external forces ``F_{x,i}``, ``F_{y,i}``, ``F_{z,i}``,
``M_{x,i}``, ``M_{y,i}``, ``M_{z,i}`` or displacements ``u_{x,i}``, ``u_{y,i}``,
``u_{z,i}``, ``\\theta_{x,i}``, ``\\theta_{y,i}``, ``\\theta_{z,i}`` applied to
each node.

** When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface **
"""
couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(aero::LiftingLine, stru::GEBT, dyn::RigidBody)
    return 5 + 6*length(stru.icol_pt)
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, dyn::RigidBody, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nd = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    id = nλ + nq + 1 : nλ + nq + nd # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    d = view(u, id) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = d
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+npd+1], p[npa+nps+npd+2], p[npa+nps+npd+3])
    ρ = p[npa+nps+npd+4]
    g = p[npa+nps+npd+5]

    # calculate gravity vector
    sϕ, cϕ = sincos(ϕr)
    sθ, cθ = sincos(θr)
    gvec = g*[-sθ, cθ*sϕ, cθ*cϕ]

    # add freestream velocity due to translation
    Vinf -= V

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start, stru.element_stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # initialize rigid body properties
    mass = 0.0
    Ir = @SMatrix zeros(3,3)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = 6*(ip-1)
        poff = npa + nps + npd + 5 + 6*(ip-1)
        ys[yoff+1] = p[poff+1]
        ys[yoff+2] = p[poff+2]
        ys[yoff+3] = p[poff+3]
        ys[yoff+4] = p[poff+4]
        ys[yoff+5] = p[poff+5]
        ys[yoff+6] = p[poff+6]
    end

    # loop through each beam element / lifting line section
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic states
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling

        # location of element (in body frame)
        pelem = element.x + u_elem

        # element length
        ΔL = element.L

        # element inertial properties
        poff = 3*npoint + 36*(i-1)
        μ = p[poff + 31]
        xm2 = p[poff + 32]
        xm3 = p[poff + 33]
        i22 = p[poff + 34]
        i33 = p[poff + 35]
        i23 = p[poff + 36]

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add local freestream linear velocities due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # calculate lifting line section inputs
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs(aero.models[i], LiftingLineSection(), ui, pi, t)

        # save aerodynamic inputs for this element
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # gravitational distributed loads (in body frame)
        fig = μ*gvec
        mig = cross(CtCab*R'*SVector(0, xm2, xm3), fig)

        # aerodynamic distributed loads (in body frame)
        fia = CtCab*R'*SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mia = CtCab*R'*SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # total distributed loads (in body frame)
        fib = fig + fia
        mib = mig + mia

        # element mass and section inertia matrix (in body frame about element center)
        melem = ΔL*μ
        Ielem = ΔL*(CtCab*R')*(@SMatrix [i22+i33 0 0; 0 i22 -i23; 0 -i23 i33])*(R*CtCab')

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        ys[yoff+1 : yoff+3] = fib
        ys[yoff+4 : yoff+6] = mib

        # add contribution to total force and moment per unit length
        Ftot += ΔL*fib
        Mtot += cross(pelem, ΔL*fib) + ΔL*mib

        # add element mass and inertia to total mass and inertia
        mass += melem
        Ir += Ielem + melem*(pelem'*pelem*I - pelem*pelem')

    end

    # save body frame linear and angular velocities
    yoff = 6*npoint + 6*nelem
    ys[yoff+1:yoff+3] = V
    ys[yoff+4:yoff+6] = Ω

    # save rigid body model inputs
    yd[1] = mass
    yd[2] = Ixx = Ir[1,1]
    yd[3] = Iyy = Ir[2,2]
    yd[4] = Izz = Ir[3,3]
    yd[5] = Ixz = Ir[1,3]
    yd[6] = Ixy = Ir[1,2]
    yd[7] = Iyz = Ir[2,3]
    yd[8:10] = Ftot
    yd[11:13] = Mtot

    return y
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT,
    dyn::RigidBody, u, p, t) where {N,T}

    # start with zero valued mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nd = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    id = nλ + nq + 1 : nλ + nq + nd # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ps = view(p, ips) # structural parameters
    d = view(u, id) # rigid body state variables
    pd = view(p, ipd) # rigid body parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = d
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])
    ρ = p[npa+nps+4]

    # add freestream velocity due to translation
    Vinf -= V
    dVinf_dV = -SMatrix{3,3}(I)

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start, stru.element_stop)


    # initialize total forces and moments
    Ftot_dV = @SMatrix zeros(3, 3)
    Ftot_dΩ = @SMatrix zeros(3, 3)
    Mtot_dV = @SMatrix zeros(3, 3)
    Mtot_dΩ = @SMatrix zeros(3, 3)

    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract aerodynamic state variables corresponding to this element
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling

        # location of element (in body frame)
        pelem = element.x + u_elem

        # element length
        ΔL = element.L

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add local freestream linear velocities due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        dv_dV = R*CtCab'*dVinf_dV
        dv_dΩ = R*CtCab'*GXBeam.tilde(pelem)
        dv_dP = -R * element.minv11 * stru.mass_scaling
        dv_dH = -R * element.minv12 * stru.mass_scaling

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        dω_dΩ = R*CtCab'
        dω_dP = R * element.minv12' * stru.mass_scaling
        dω_dH = R * element.minv22 * stru.mass_scaling

        # calculate lifting line section input mass matrix
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        Myi = get_input_mass_matrix(aero.models[i], LiftingLineSection(), ui, pi, t)

        # separate into component mass matrices
        d_dλ = SMatrix{Nyai[i], Nλi[i]}(view(Myi, 1:Nyai[i], 1:Nλi[i]))

        d_dv = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+1 : Nλi[i]+3))
        d_dω = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+4 : Nλi[i]+6))
        d_dV = d_dv * dv_dV
        d_dΩ = d_dv * dv_dΩ + d_dω * dω_dΩ
        d_dP = d_dv * dv_dP + d_dω * dω_dP
        d_dH = d_dv * dv_dH + d_dω * dω_dH

        f_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+1 : Nyai[i]+3, 1:Nλi[i]))
        m_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+4 : Nyai[i]+6, 1:Nλi[i]))

        f_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+1 : Nλi[i]+3))
        f_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+4 : Nλi[i]+6))
        f_dV = f_dv * dv_dV
        f_dΩ = f_dv * dv_dΩ + f_dω * dω_dΩ
        f_dP = f_dv * dv_dP + f_dω * dω_dP
        f_dH = f_dv * dv_dH + f_dω * dω_dH

        m_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+1 : Nλi[i]+3))
        m_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+4 : Nλi[i]+6))
        m_dV = m_dv * dv_dV
        m_dΩ = m_dv * dv_dΩ + m_dω * dω_dΩ
        m_dP = m_dv * dv_dP + m_dω * dω_dP
        m_dH = m_dv * dv_dH + m_dω * dω_dH

        # save aerodynamic input mass matrix entries (in local aerodynamic frame)
        icol = stru.icol_beam[i]
        My[iyas[i], iλs[i]] = d_dλ
        My[iyas[i], nλ+icol+12:nλ+icol+14] = d_dP
        My[iyas[i], nλ+icol+15:nλ+icol+17] = d_dH
        My[iyas[i], nλ+nq+7 : nλ+nq+9] = d_dV
        My[iyas[i], nλ+nq+10 : nλ+nq+12] = d_dΩ

        # extract local forces and moments (in body frame)
        fi_dλ = CtCab*R'*f_dλ
        fi_dP = CtCab*R'*f_dP
        fi_dH = CtCab*R'*f_dH
        fi_dV = CtCab*R'*f_dV
        fi_dΩ = CtCab*R'*f_dΩ
        mi_dλ = CtCab*R'*m_dλ
        mi_dP = CtCab*R'*m_dP
        mi_dH = CtCab*R'*m_dH
        mi_dV = CtCab*R'*m_dV
        mi_dΩ = CtCab*R'*m_dΩ

        # save load mass matrix entries (in body frame)
        offset = nya + 6*(i-1)
        My[offset+1 : offset+3, iλs[i]] = fi_dλ
        My[offset+4 : offset+6, iλs[i]] = mi_dλ
        My[offset+1 : offset+3, nλ+icol+12:nλ+icol+14] = fi_dP
        My[offset+1 : offset+3, nλ+icol+15:nλ+icol+17] = fi_dH
        My[offset+4 : offset+6, nλ+icol+12:nλ+icol+14] = mi_dP
        My[offset+4 : offset+6, nλ+icol+15:nλ+icol+17] = mi_dH
        My[offset+1 : offset+3, nλ+nq+7 : nλ+nq+9] = fi_dV
        My[offset+1 : offset+3, nλ+nq+10 : nλ+nq+12] = fi_dΩ
        My[offset+4 : offset+6, nλ+nq+7 : nλ+nq+9] = mi_dV
        My[offset+4 : offset+6, nλ+nq+10 : nλ+nq+12] = mi_dΩ

        # add contribution to total forces and moments
        My[nya+nys+8 : nya+nys+10, iλs[i]] = ΔL*fi_dλ
        My[nya+nys+8 : nya+nys+10, nλ+icol+12:nλ+icol+14] = ΔL*f_dP
        My[nya+nys+8 : nya+nys+10, nλ+icol+15:nλ+icol+17] = ΔL*f_dH
        Ftot_dV += ΔL*fi_dV
        Ftot_dΩ += ΔL*fi_dΩ

        My[nya+nys+11 : nya+nys+13, iλs[i]] = GXBeam.tilde(pelem)*ΔL*fi_dλ + ΔL*mi_dλ
        My[nya+nys+11 : nya+nys+13, nλ+icol+12:nλ+icol+14] = ΔL*mi_dP
        My[nya+nys+11 : nya+nys+13, nλ+icol+15:nλ+icol+17] = ΔL*mi_dH
        Mtot_dV += GXBeam.tilde(pelem)*ΔL*fi_dV + ΔL*mi_dV
        Mtot_dΩ += GXBeam.tilde(pelem)*ΔL*fi_dΩ + ΔL*mi_dΩ

    end

    # save rigid body model inputs
    My[nya+nys+8 : nya+nys+10, nλ+nq+7:nλ+nq+9] = Ftot_dV
    My[nya+nys+8 : nya+nys+10, nλ+nq+10:nλ+nq+12] = Ftot_dΩ
    My[nya+nys+11 : nya+nys+13, nλ+nq+7:nλ+nq+9] = Mtot_dV
    My[nya+nys+11 : nya+nys+13, nλ+nq+10:nλ+nq+12] = Mtot_dΩ

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::GEBT, dyn::RigidBody, du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Nu = number_of_inputs(models)
    y = zeros(TF, Nu)

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nd = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    id = nλ + nq + 1 : nλ + nq + nd # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dλ = view(du, iλ) # aerodynamic states rates
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dq = view(du, iq) # structural state rates
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    dd = view(du, id) # rigid body state rates
    d = view(u, id) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    dλs = view.(Ref(dλ), iλs) # aerodynamic state rates for each section
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = d
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # extract freestream parameters
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])
    ρ = p[npa+nps+4]

    # add freestream velocity due to translation
    Vinf -= V

    # number of points and elements
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.element_start, stru.element_stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic state rates and states
        dλi = SVector{Nλi[i]}(dλs[i])
        λi = SVector{Nλi[i]}(λs[i])

        # get element structural state rates and states
        icol = stru.icol_beam[i]
        u_elem = SVector(q[icol], q[icol+1], q[icol+2])
        θ_elem = SVector(q[icol+3], q[icol+4], q[icol+5])
        F_elem = SVector(q[icol+6], q[icol+7], q[icol+8]) .* stru.force_scaling
        M_elem = SVector(q[icol+9], q[icol+10], q[icol+11]) .* stru.force_scaling
        P_elem = SVector(q[icol+12], q[icol+13], q[icol+14]) .* stru.mass_scaling
        H_elem = SVector(q[icol+15], q[icol+16], q[icol+17]) .* stru.mass_scaling

        du_elem = SVector(dq[icol], dq[icol+1], dq[icol+2])
        dθ_elem = SVector(dq[icol+3], dq[icol+4], dq[icol+5])
        dF_elem = SVector(dq[icol+6], dq[icol+7], dq[icol+8]) .* stru.force_scaling
        dM_elem = SVector(dq[icol+9], dq[icol+10], dq[icol+11]) .* stru.force_scaling
        dP_elem = SVector(dq[icol+12], dq[icol+13], dq[icol+14]) .* stru.mass_scaling
        dH_elem = SVector(dq[icol+15], dq[icol+16], dq[icol+17]) .* stru.mass_scaling

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling
        dθ_elem *= scaling

        # element length
        ΔL = element.L

        # location of element (in body frame)
        pelem = element.x + u_elem

        # transformation matrix from local beam frame to the body frame

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        Ct = GXBeam.get_C(θ_elem)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local beam frame to local aerodynamic frame

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from body frame to local aerodynamic frame
        vi = R*CtCab'*Vinf

        # add velocity due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, P_elem, H_elem)

        # calculate freestream acceleration due to surface motion
        dvi = -R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)

        # calculate local section angular velocities due to surface motion
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # calculate local section angular accelerations due to surface motion
        dωi = R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # calculate lifting line section inputs from state rate contributions
        dui = vcat(dλi, dvi, dωi) # section state rates
        ui = vcat(λi, vi, ωi) # section states
        pi = vcat(SVector{Npai[i]}(pas[i]), ρ) # section parameters
        yi = get_inputs_from_state_rates(aero.models[i], LiftingLineSection(),
            dui, ui, pi, t)

        # extract and save aerodynamic model inputs (in local aerodynamic frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments per unit length (in local aerodynamic frame)
        fi = SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mi = SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # force and moment per unit length (in body frame)
        fbi = CtCab*R'*fi
        mbi = CtCab*R'*mi

        # add apparent forces due to linear and angular acceleration
        fbi += melem*dV + cross(pelem, melem*dΩ)
        mbi += melem*dΩ

        # save element force and moment per unit length (in body frame)
        ys[6*(i-1)+1 : 6*(i-1)+3] = fbi
        ys[6*(i-1)+4 : 6*(i-1)+6] = mbi

        # add contribution to total force and moment per unit length
        Ftot += ΔL*fbi
        Mtot += cross(pelem, ΔL*fbi) + ΔL*mbi

    end

    # save state rate contribution to body frame linear and angular velocities
    ys[6*N+1 : 6*N+3] .= 0
    ys[6*N+4 : 6*N+6] .= 0

    # save state rate contribution to rigid body model inputs
    yd[1] = 0
    yd[2] = Ixx = 0
    yd[3] = Iyy = 0
    yd[4] = Izz = 0
    yd[5] = Ixz = 0
    yd[6] = Ixy = 0
    yd[7] = Iyz = 0
    yd[8:10] = Ftot
    yd[11:13] = Mtot

    return y
end

# --- Internal Methods --- #

function input_jacobian_product!(y, x, model::LiftingLine)

    models = model.models

    iy = state_indices(models)
    ix = input_indices(models)

    Jyi = get_input_jacobian.(models)

    yi = view.(Ref(y), iy)
    xi = view.(Ref(x), ix)

    mul!.(yi, Jyi, xi)

    return y
end

function input_jacobian_product!(y, x, model::LiftingLine, λ, d, p, t)

    models = model.models

    iu = state_indices(models)
    id = input_indices(models)
    ip = parameter_indices(models)

    xi = view.(Ref(x), id)
    yi = view.(Ref(y), iu)

    λi = view.(Ref(λ), iu)
    di = view.(Ref(d), id)
    pi = view.(Ref(p), ip)

    Jyi = get_input_jacobian.(models, λi, di, pi, t)

    mul!.(yi, Jyi, xi)

    return y
end
