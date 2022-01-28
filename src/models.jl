"""
    Submodel{iip, F, NS, NI, NP, JR, JS, JI, JP, GT, IS, II, IP, OS, OI, OP}

Base type for submodels in `Aeroelasticity.jl`.
"""
struct Submodel{iip, F, NS, NI, NP, C, JR, JS, JI, JP, GT, IS, II, IP, OS, OI, OP}
    # residual function
    f::F
    # dimensions
    nx::NS
    ny::NI
    np::NP
    # constants
    constants::C
    # jacobians/gradients
    ratejac::JR
    statejac::JS
    inputjac::JI
    paramjac::JP
    tgrad::GT
    # input/output functions
    setstate::IS
    setinput::II
    setparam::IP
    sepstate::OS
    sepinput::OI
    sepparam::OP
end

"""
    Submodel{iip}(f, nx, ny, np; kwargs...)

Define a custom submodel for use with `Aeroelasticity.jl`.  Note that `iip` should be set to
`true` if `f` is defined in-place and `false` otherwise. 

# Arguments
 - `f`: Residual function (`resid = f(dx, x, y, p, t)` or `f(resid, dx, x, y, p, t)`)
 - `nx`: Number of states.
 - `ny`: Number of inputs.
 - `np`: Number of parameters.

# Keyword Arguments
 - `ratejac = Nonlinear()`:  Rate jacobian description (of type [`AbstractJacobian`](@ref))
 - `statejac = Nonlinear()`: State jacobian description (of type [`AbstractJacobian`](@ref))
 - `inputjac = Nonlinear()`: Input jacobian description (of type [`AbstractJacobian`](@ref))
 - `paramjac = Nonlinear()`: Parameter jacobian description (of type [`AbstractJacobian`](@ref))
 - `tgrad = Nonlinear()`: Time gradient description (of type [`AbstractGradient`](@ref))
 - `setstate = (states; x) -> states .= x`: Function for setting rates/states using a 
    custom set of keyword arguments.
 - `setinput = (inputs; y) -> inputs .= y`: Function for setting inputs using a custom set 
    of keyword arguments.
 - `setparam = (params; p) -> params .= p`: Function for setting parameters using a custom 
    set of keyword arguments.
 - `sepstate = (x) -> (x = x,)`: Function for separating rates/states into named variables.
 - `sepinput = (y) -> (y = y,)`: Function for separating inputs into named variables.
 - `sepparam = (p) -> (p = p,)`: Function for separating parameters into named variables.
"""
function Submodel{iip}(f, nx, ny, np;
    constants = (),
    ratejac = Nonlinear(),
    statejac = Nonlinear(),
    inputjac = Nonlinear(),
    paramjac = Nonlinear(),
    tgrad = Nonlinear(),
    setstate = (states; x) -> states .= x,
    setinput = (inputs; y) -> inputs .= y,
    setparam = (params; p) -> params .= p,
    sepstate = (x) -> (x = x,),
    sepinput = (y) -> (y = y,),
    sepparam = (p) -> (p = p,)) where iip

    # add jacobian definitions, if necessary
    ratejac = add_rate_jacobian(ratejac, f, iip, nx, ny, np)
    statejac = add_state_jacobian(statejac, f, iip, nx, ny, np)
    inputjac = add_input_jacobian(inputjac, f, iip, nx, ny, np)
    paramjac = add_parameter_jacobian(paramjac, f, iip, nx, ny, np)
    tgrad = add_time_gradient(tgrad, f, iip, nx, ny, np)

    return Submodel{iip,
        typeof(f), typeof(nx), typeof(ny), typeof(np), typeof(constants),
        typeof(ratejac), typeof(statejac), typeof(inputjac), typeof(paramjac), typeof(tgrad),
        typeof(setstate), typeof(setinput), typeof(setparam),
        typeof(sepstate), typeof(sepinput), typeof(sepparam)}(
            f, nx, ny, np, constants,
            ratejac, statejac, inputjac, paramjac, tgrad,
            setstate, setinput, setparam, 
            sepstate, sepinput, sepparam)
end

"""
    Submodel(np; kwargs...)

Define a submodel with no state variables (see [`Submodel`](@ref))

# Arguments
 - `np`: Number of parameters.

# Keyword Arguments
 - `constants = ()`: Constants associated with this model
 - `setparam = (params; p) -> params .= p`: Function for setting parameters using a custom 
    set of keyword arguments.
 - `sepparam = (p) -> (p = p,)`: Function for separating parameters into named variables.
"""
function Submodel(np; 
    constants = (),
    setparam = (params; p) -> params .= p,
    sepparam = (p) -> (p = p,))

    return Submodel{false}((dx, x, y, p, t) -> Float64[], Val(0), Val(0), np;
        constants = constants,
        ratejac = Empty(),
        statejac = Empty(),
        inputjac = Empty(),
        paramjac = Empty(),
        tgrad = Empty(),
        setstate = (x) -> x,
        setinput = (y) -> y,
        setparam = setparam,
        sepstate = (x) -> (),
        sepinput = (y) -> (),
        sepparam = sepparam)
end

"""
    Coupling{iip, G, NX, NY, NP, NPC, JR, JS, JP, GT, IP, OP}

Base type for coupling models in Aeroelasticity.jl.
"""
struct Coupling{iip, G, NX, NY, NP, NPC, JR, JS, JP, GT, IP, OP}
    # coupling function
    g::G
    # dimensions
    nx::NX
    ny::NY
    np::NP
    npc::NPC
    # jacobians/gradients
    ratejac::JR
    statejac::JS
    paramjac::JP
    tgrad::GT
    # input/output functions
    setparam::IP
    sepparam::OP
end

"""
    Coupling{iip}(g, nx, ny, np, npc; kwargs...)

Define a custom coupling model for use with Aeroelasticity.jl.  Note that `iip` should 
be set to `true` if `g` is defined in-place and `false` otherwise. 

# Arguments
 - `g`: Coupling function (`y = g(dx, x, p, t)` or `g(y, dx, x, p, t)`)
 - `nx`: Number of states.
 - `ny`: Number of inputs.
 - `np`: Number of parameters (including parameters introduced by the coupling).
 - `npc`: Number of additional parameters introduced by the coupling.

# Keyword Arguments
 - `ratejac = Nonlinear()`:  Rate jacobian description (of type [`AbstractJacobian`](@ref))
 - `statejac = Nonlinear()`: State jacobian description (of type [`AbstractJacobian`](@ref))
 - `paramjac = Nonlinear()`: Parameter jacobian description (of type [`AbstractJacobian`](@ref))
 - `tgrad = Nonlinear()`: Time gradient description (of type [`AbstractGradient`](@ref))
 - `setparam = (params; p) -> params .= p`: Function for setting parameters introduced by
    the coupling using a custom set of keyword arguments.
 - `sepparam = (p) -> (p = p,)`: Function for separating parameters introduced by the 
    coupling into named variables.
"""
function Coupling{iip}(g, nx, ny, np, npc;
    ratejac = Nonlinear(),
    statejac = Nonlinear(),
    paramjac = Nonlinear(),
    tgrad = Nonlinear(),
    setparam = (params; p) -> params .= p,
    sepparam = (p) -> (p = p,)) where {iip}

    # add jacobian definitions, if necessary
    ratejac = add_coupling_rate_jacobian(ratejac, g, iip, nx, ny, np)
    statejac = add_coupling_state_jacobian(statejac, g, iip, nx, ny, np)
    paramjac = add_coupling_parameter_jacobian(paramjac, g, iip, nx, ny, np)
    tgrad = add_coupling_time_gradient(tgrad, g, iip, nx, ny, np)

    return Coupling{iip, 
        typeof(g), typeof(nx), typeof(ny), typeof(np), typeof(npc),
        typeof(ratejac), typeof(statejac), typeof(paramjac), typeof(tgrad),
        typeof(setparam), typeof(sepparam)}(
            g, nx, ny, np, npc, 
            ratejac, statejac, paramjac, tgrad,
            setparam, sepparam)
end

"""
    CoupledModel{M,C}

Base type for coupled models in Aeroelasticity.jl.
"""
struct CoupledModel{M,C}
    submodels::M
    coupling::C
end

"""
    CoupledModel(submodels, coupling)

Define a coupled model given a tuple of submodels (of type [`Submodel`](@ref)) and a 
coupling model (of type [`Coupling`](@ref)).
"""
CoupledModel(submodels, coupling)

# --- Inplaceness Trait --- #

abstract type InPlaceness end
struct InPlace <: InPlaceness end
struct OutOfPlace <: InPlaceness end

inplaceness(::Submodel{true}) = InPlace()
inplaceness(::Submodel{false}) = OutOfPlace()

inplaceness(::Coupling{true}) = InPlace()
inplaceness(::Coupling{false}) = OutOfPlace()

function inplaceness(model::CoupledModel)
    if any(isinplace.(model.submodels)) || isinplace(model.coupling)
        return InPlace()
    else
        return OutOfPlace()
    end
end

isinplace(x) = isinplace(inplaceness(x))
isinplace(::InPlace) = true
isinplace(::OutOfPlace) = false