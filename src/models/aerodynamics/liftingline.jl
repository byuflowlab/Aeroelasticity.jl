"""
    LiftingLine

Lifting line model constructed by concatenating the governing equations, state variables,
inputs, and parameters of multiple two-dimensional aerodynamic models.
"""
struct LiftingLine{T}
    section_models::T
end

"""
    LiftingLine(section_models)

Construct a lifting line aerodynamic model from a collection of two-dimensional aerodynamic 
models.
"""
LiftingLine(section_models)

# --- Submodel Creation --- #

function Submodel(model::LiftingLine)

    section_models = model.section_models

    # residual function
    fresid = (resid, dx, x, y, p, t) -> liftingline_residual!(resid, dx, x, y, p, t; 
        section_models = section_models)

    # number of states, inputs, and parameters
    nx = sum(number_of_states.(section_models))
    ny = sum(number_of_inputs.(section_models))
    np = sum(number_of_parameters.(section_models))

    # rate jacobian definition (depends on submodel rate jacobian definitions)
    rate_jacobians = getproperty.(section_models, :ratejac)
    if all(isempty.(rate_jacobians))
        ratejac = Empty()
    elseif all(iszero.(rate_jacobians))
        ratejac = Zeros()
    elseif all(isidentity.(rate_jacobians))
        ratejac = Identity()
    elseif all(isinvariant.(rate_jacobians))
        ratejac = Invariant(liftingline_rate_jacobian(; section_models))
    elseif all(isconstant.(rate_jacobians))
        ratejac = Constant((J, p) -> liftingline_rate_jacobian!(J, p; section_models))
    elseif all(islinear.(rate_jacobians))
        ratejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; section_models)
        )
    else
        ratejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; section_models)
        )
    end

    # state jacobian definition (depends on submodel state jacobian definitions)
    state_jacobians = getproperty.(section_models, :statejac)
    if all(isempty.(state_jacobians))
        statejac = Empty()
    elseif all(iszero.(state_jacobians))
        statejac = Zeros()
    elseif all(isidentity.(state_jacobians))
        statejac = Identity()
    elseif all(isinvariant.(state_jacobians))
        statejac = Invariant(liftingline_state_jacobian(; section_models))
    elseif all(isconstant.(state_jacobians))
        statejac = Constant((J, p) -> liftingline_state_jacobian!(J, p; section_models))
    elseif all(islinear.(state_jacobians))
        statejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; section_models)
        )
    else
        statejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; section_models)
        )
    end

    # input jacobian definition (depends on submodel input jacobian definitions)
    input_jacobians = getproperty.(section_models, :inputjac)
    if all(isempty.(input_jacobians))
        inputjac = Empty()
    elseif all(iszero.(input_jacobians))
        inputjac = Zeros()
    elseif all(isidentity.(input_jacobians))
        inputjac = Identity()
    elseif all(isinvariant.(input_jacobians))
        inputjac = Invariant(liftingline_input_jacobian(; section_models))
    elseif all(isconstant.(input_jacobians))
        inputjac = Constant((J, p) -> liftingline_input_jacobian!(J, p; section_models))
    elseif all(islinear.(input_jacobians))
        inputjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; section_models)
        )
    else
        inputjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; section_models)
        )
    end

    # parameter jacobian definition (depends on submodel parameter jacobian definitions)
    parameter_jacobians = getproperty.(section_models, :paramjac)
    if all(isempty.(parameter_jacobians))
        paramjac = Empty()
    elseif all(iszero.(parameter_jacobians))
        paramjac = Zeros()
    elseif all(isidentity.(parameter_jacobians))
        paramjac = Identity()
    elseif all(isinvariant.(parameter_jacobians))
        paramjac = Invariant(liftingline_parameter_jacobian(; section_models))
    elseif all(isconstant.(parameter_jacobians))
        paramjac = Constant((J, p) -> liftingline_parameter_jacobian!(J, p; section_models))
    elseif all(islinear.(parameter_jacobians))
        paramjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; section_models)
        )
    else
        paramjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; section_models)
        )
    end

    # time gradient definition (depends on submodel time gradient definitions)
    time_gradients = getproperty.(section_models, :tgrad)
    if all(isempty.(time_gradients))
        tgrad = Empty()
    elseif all(iszero.(time_gradients))
        tgrad = Zeros()
    elseif all(isinvariant.(time_gradients))
        tgrad = Invariant(liftingline_time_gradient(; section_models))
    elseif all(isconstant.(time_gradients))
        tgrad = Constant((dT, p) -> liftingline_time_gradient!(dT, p; section_models))
    elseif all(islinear.(time_gradients))
        paramjac = Linear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; section_models)
        )
    else
        paramjac = Nonlinear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; section_models)
        )
    end
    
    # convenience functions for setting states, inputs, and parameters
    setstate = (x; kwargs...) -> liftingline_setstate!(x; section_models, kwargs...)
    setinput = (y; kwargs...) -> liftingline_setinput!(y; section_models, kwargs...)
    setparam = (p; kwargs...) -> liftingline_setparam!(p; section_models, kwargs...)
    
    # convenience functions for separating states, inputs, and parameters
    sepstate = (x) -> liftingline_sepstate(x; section_models)
    sepinput = (y) -> liftingline_sepinput(y; section_models)
    sepparam = (p) -> liftingline_sepparam(p; section_models)

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
function liftingline_residual!(resid, dx, x, y, p, t; section_models)

    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)

    resids = view.(Ref(resid), ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_residual!.(resids, section_models, dxs, xs, ys, ps, t)

    return resid
end

# rate jacobian function (for an invariant rate jacobian)
function liftingline_rate_jacobian(; section_models)
    
    nx = sum(number_of_states.(section_models))

    J = zeros(nx, nx)
    
    ix = state_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    
    get_rate_jacobian!.(Js, section_models)
    
    return J
end

# rate jacobian function (for a constant rate jacobian)
function liftingline_rate_jacobian!(J, p; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_rate_jacobian!.(Js, section_models, ps)
    
    return J
end

# rate jacobian function (for a time-dependent rate jacobian)
function liftingline_rate_jacobian!(J, dx, x, y, p, t; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_rate_jacobian!.(Js, section_models, dxs, xs, ys, ps, t)
    
    return J
end

# state jacobian function (for an invariant state jacobian)
function liftingline_state_jacobian(; section_models)

    nx = sum(number_of_states.(section_models))

    J = zeros(nx, nx)
    
    ix = state_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    
    get_state_jacobian!.(Js, section_models)

    return J
end

# state jacobian function (for a constant state jacobian)
function liftingline_state_jacobian!(J, p; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_state_jacobian!.(Js, section_models, ps)
    
    return J
end

# state jacobian function (for a time-dependent state jacobian)
function liftingline_state_jacobian!(J, dx, x, y, p, t; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_state_jacobian!.(Js, section_models, dxs, xs, ys, ps, t)
    
    return J
end

# input jacobian function (for an invariant input jacobian)
function liftingline_input_jacobian(; section_models)

    nx = sum(number_of_states.(section_models))
    ny = sum(number_of_inputs.(section_models))

    J = zeros(nx, ny)
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    
    Js = view.(Ref(J), ix, iy)
    
    get_input_jacobian!.(Js, section_models)

    return J
end

# input jacobian function (for a constant input jacobian)
function liftingline_input_jacobian!(J, p; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, iy)
    ps = view.(Ref(p), ip)
    
    get_input_jacobian!.(Js, section_models, ps)
    
    return J
end

# input jacobian function (for a time-dependent input jacobian)
function liftingline_input_jacobian!(J, dx, x, y, p, t; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, iy)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_input_jacobian!.(Js, section_models, dxs, xs, ys, ps, t)
    
    return J
end

# parameter jacobian function (for an invariant parameter jacobian)
function liftingline_parameter_jacobian(; section_models)
    
    nx = sum(number_of_states.(section_models))
    np = sum(number_of_parameters.(section_models))

    J = zeros(nx, np)
    
    ix = state_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ip)
    
    get_parameter_jacobian!.(Js, section_models)

    return J
end

# parameter jacobian function (for a constant parameter jacobian)
function liftingline_parameter_jacobian!(J, p; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ip)
    ps = view.(Ref(p), ip)
    
    get_parameter_jacobian!.(Js, section_models, ps)
    
    return J
end

# parameter jacobian function (for a time-dependent parameter jacobian)
function liftingline_parameter_jacobian!(J, dx, x, y, p, t; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    Js = view.(Ref(J), ix, ip)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_parameter_jacobian!.(Js, section_models, dxs, xs, ys, ps, t)
    
    return J
end

# time gradient function (for an invariant time gradient)
liftingline_time_gradient(; section_models) = vcat(get_time_gradient.(section_models)...)

# time gradient function (for a constant time gradient)
function liftingline_time_gradient!(dT, p; section_models)

    J .= 0
    
    ix = state_indices(section_models)
    ip = parameter_indices(section_models)
    
    dTs = view.(Ref(dT), ix)
    ps = view.(Ref(p), ip)
    
    get_time_gradient!.(dTs, section_models, ps)
    
    return dT
end

# time gradient function (for a time-dependent time gradient)
function liftingline_time_gradient!(dT, dx, x, y, p, t; section_models)

    dT .= 0
    
    ix = state_indices(section_models)
    iy = input_indices(section_models)
    ip = parameter_indices(section_models)
    
    dTs = view.(Ref(dT), ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_time_gradient!.(dTs, section_models, dxs, xs, ys, ps, t)
    
    return J
end

# convenience function for defining this model's state vector
function liftingline_setstate!(x; sections, section_models)
    section_indices = state_indices(section_models)
    vxs = view.(Ref(x), section_indices)
    bfn! = (x, model, states) -> set_states!(x, model; states...)
    bfn!.(vxs, section_models, sections)    
    return x
end

# convenience function for defining this model's input vector
function liftingline_setinput!(y; sections, section_models)
    section_indices = input_indices(section_models)
    vys = view.(Ref(y), section_indices)
    bfn! = (y, model, inputs) -> set_inputs!(y, model; inputs...)
    bfn!.(vys, section_models, sections)
    return y
end

# convenience function for defining this model's parameter vector
function liftingline_setparam!(p; sections, section_models)
    section_indices = parameter_indices(section_models)
    vps = view.(Ref(p), section_indices)
    bfn! = (p, model, parameters) -> set_parameters!(p, model; parameters...)
    bfn!.(vps, section_models, sections)
    return p
end

# convenience function for separating this model's state vector
function liftingline_sepstate(x; section_models)
    section_indices = state_indices(section_models)
    vxs = view.(Ref(x), section_indices)
    return (sections = separate_states.(section_models, vxs),)
end

# convenience function for separating this model's input vector
function liftingline_sepinput(y; section_models) 
    section_indices = input_indices(section_models)
    vys = view.(Ref(y), section_indices)
    return (sections = separate_inputs.(section_models, vys),)
end

# convenience function for separating this model's parameter vector
function liftingline_sepparam(p; section_models)
    section_indices = parameter_indices(section_models)
    vps = view.(Ref(p), section_indices)
    return (sections = separate_parameters.(section_models, vps),)
end

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