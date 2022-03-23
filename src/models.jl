"""
    Submodel{iip, NS, NI, NP}

Base type for submodels in `Aeroelasticity.jl`.
"""
struct Submodel{iip, NS, NI, NP}
    # residual function
    f
    # dimensions
    nx::NS
    ny::NI
    np::NP
    # jacobians/gradients
    ratejac
    statejac
    inputjac
    paramjac
    tgrad
    # input/output functions
    setstate
    setinput
    setparam
    sepstate
    sepinput
    sepparam
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

    return Submodel{iip, typeof(nx), typeof(ny), typeof(np)}(
            f, nx, ny, np, 
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
 - `setparam = (params; p) -> params .= p`: Function for setting parameters using a custom 
    set of keyword arguments.
 - `sepparam = (p) -> (p = p,)`: Function for separating parameters into named variables.
"""
function Submodel(np; 
    setparam = (params; p) -> params .= p,
    sepparam = (p) -> (p = p,))

    return Submodel{false}((dx, x, y, p, t) -> Float64[], Val(0), Val(0), np;
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
    Coupling{iip, G, NX, NY, NP, NPC, JR, JS, JP, GT}

Base type for coupling models in Aeroelasticity.jl.
"""
struct Coupling{iip, NX, NY, NP, NPC}
    # coupling function
    g
    # dimensions
    nx::NX
    ny::NY
    np::NP
    npc::NPC
    # jacobians/gradients
    ratejac
    statejac
    paramjac
    tgrad
    # input/output functions
    setparam
    sepparam
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

    return Coupling{iip, typeof(nx), typeof(ny), typeof(np), typeof(npc)}(
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