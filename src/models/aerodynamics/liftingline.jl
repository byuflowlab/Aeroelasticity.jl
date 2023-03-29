"""
    LiftingLine{R, Y, I}

Lifting line model constructed by concatenating the governing equations, state variables,
inputs, and parameters of multiple two-dimensional aerodynamic models.
"""
struct LiftingLine{R, Y, I}
    fresid::R
    finput::Y
    indices::I
end

"""
    LiftingLine(fresid, nstate=number_of_states.(fresid); finput=default_coupling.(fresid, LiftingLineSection()))

Initializes a lifting line model of type [`LiftingLine`](@ref)

 # Arguments:
  - `fresid`: Residual function for each lifting line section
  - `nstate`: Number of state variables for each lifting line section

 # Keyword Arguments:
  - `finput`: Coupling function for each lifting line section.
"""
function LiftingLine(fresid, nstate=number_of_states.(fresid);
    finput=default_coupling.(fresid, Ref(LiftingLineSection())))

    # compute state indices
    ix2 = cumsum(nstate)
    ix1 = ix2 .- nstate .+ 1
    indices = UnitRange.(ix1, ix2)

    # determine types
    TR = typeof(fresid)
    TY = typeof(finput)
    TI = typeof(indices)

    return LiftingLine{TR, TY, TI}(fresid, finput, indices)
end

# residual function
function (liftingline::LiftingLine)(resid, dx, x, y, p, t)

    # extract constants
    @unpack fresid, indices = liftingline

    # separate residuals
    residuals = view.(Ref(resid), indices)

    # separate rates
    rates = view.(Ref(dx), indices)

    # separate states
    states = view.(Ref(x), indices)

    # inputs are already separated
    inputs = y

    # parameters are already separated
    parameters = p

    # compute residuals (in-place)
    for i in eachindex(fresid)
        fresid[i](residuals[i], rates[i], states[i], inputs[i], parameters[i], t)
    end

    # return modified residual vector
    return resid
end

number_of_states(liftingline::LiftingLine) = liftingline.indices[end][end]

"""
    LiftingLineParameters(section_parameters)

Defines parameters for a lifting line theory aerodynamic model.
"""
LiftingLineParameters(section_parameters) = section_parameters

# --- Internal Methods for Couplings with this Model --- #

# local section velocities/accelerations
function liftingline_section_velocities(x, e1, e2, e3, V, Ω, a, α)

    # transformation from body to local frame
    R = [e1 e2 e3]'

    # linear and angular freestream velocity
    vi = -R*(V + cross(Ω, x))
    ωi = R*Ω

    # linear and angular freestream acceleration
    ai = -R*(a + cross(α, x))
    αi = R*α

    return vi, ωi, ai, αi
end

function liftingline_section_inputs(model, dxa, xa, pa, v, ω, a, α, rho, c, t)

    # aerodynamic rates, states, and parameters
    dxa = dxa
    xa = xa
    pa = pa

    # structural rates, states, and parameters
    dxs = SVector(a[1], a[2], a[3], α[1], α[2], α[3])
    xs = SVector(v[1], v[2], v[3], ω[1], ω[2], ω[3])
    ps = SVector(rho, c)

    # combined rates, states, and parameters
    dx = vcat(dxa, dxs)
    x = vcat(xa, xs)
    p = vcat(pa, ps)

    # calculate section inputs
    return get_coupling_inputs(model, dx, x, p, t)
end