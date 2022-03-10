# Developer Guide

This package has been designed to extensible so that custom models may be implemented and used as necessary.  This guide details the process for defining and coupling new models.  

## Defining a New Submodel

New submodels may be defined using the following constructor:

```@docs
Aeroelasticity.Submodel{iip}(f, nx, ny, np; kwargs...)
```

The first argument, `f`, is a residual function which may be used to define the submodel's state rates ``\dot{x}`` as a function of the model's states ``x``, inputs ``y``, and parameters ``p``, as well as the current time ``t``.  Note that for the purposes of this package we consider all time varying parameters to be inputs, rather than parameters.  The function signature for this function is either 
`f!(resid, dx, x, y, p, t)` or `resid = f(dx, x, y, p, t)` depending on whether the submodel is specified as in-place (`iip=true`) or out-of-place (`iip=false`).  The arguments `dx`, `x`, `y`, and `p` are assumed to be vectors.  

The next three arguments (`nx`, `ny`, and `np`) define the number of state variables, inputs, and parameters respectively.

The remaining arguments are optional, but may be used to specify additional information about a model which can make the model easier to work with.

The gradients/jacobians associated with a model may be specified using the `ratejac`, `statejac`, `inputjac`, `paramjac`, and `tgrad` keyword arguments.  Each of these arguments, if provided, must be wrapped as an object of type [`AbstractJacobian`](@ref).

The arguments `setstate`, `setinput`, and `setparam` define convenience functions for setting the values of the state, input, and parameter vectors using keyword arguments.

The arguments `sepstate`, `sepinput`, and `sepparam` define convenience functions for separating the values of the state, input, and parameter vectors into a named tuple.

## Defining a New Model Coupling

New couplings may be defined using the following constructor:

```@docs
Aeroelasticity.Coupling{iip}(g, nx, ny, np, npc; kwargs...)
```

The first argument, `g`, is a coupling function which may be used to define the coupled model's inputs ``y`` as a function of the coupled model's state rates ``dx``, states ``x`` and parameters ``p``, as well as the current time ``t``. The function signature for this function is either `g!(y, dx, x, y, p, t)` or `y = g(dx, x, p, t)` depending on whether the model is specified as in-place (`iip=true`) or out-of-place (`iip=false`).  The arguments `dx`, `x`, and `y`, are, and `p` are assumed to be vectors.  

The next three arguments (`nx`, `ny`, and `np`) define the number of state variables, inputs, and parameters of the coupled model.  The argument `npc` defines the number of additional parameters introduced by the coupling.

The remaining arguments are optional, but may be used to specify additional information about the coupling which can make the coupled model easier to work with.

The gradients/jacobians associated with a coupling may be specified using the `ratejac`, `statejac`, `paramjac`, and `tgrad` keyword arguments.  Each of these arguments, if provided, must be wrapped as an object of type [`AbstractJacobian`](@ref).

The argument `setparam` may be provided to define a convenience function for setting the values of the parameter vector associated with the coupling function using keyword arguments.

The argument `sepparam` may be provided to define a convenience function for separating the values of the parameter vector associated with the coupling function into a named tuple.

## Defining Jacobians

The following types are available for defining gradients and jacobians:

```@docs
Aeroelasticity.AbstractJacobian
Aeroelasticity.Empty
Aeroelasticity.Zeros
Aeroelasticity.Identity
Aeroelasticity.Invariant
Aeroelasticity.Constant
Aeroelasticity.Linear
Aeroelasticity.Nonlinear
```

In general, the most specific jacobian specification should be used whenever possible, though this is most important when defining `ratejac`.  

!!! tip "Explicit Models"
    Residual functions for explicit models, (i.e. models whose governing equations can be written as ``\dot{x} = f(x, y, p, t)``) should be written in the form ``0 = \dot{x} - f(x, y, p, t)`` when used with this package.  Setting `ratejac = Identity()` will then allow this package to recover the explicit relationship from the specified residual function.

## Defining a New Coupled Model

Custom coupled models may be created by using the [`assemble_model`](@ref) function.  This function uses the following order when concatenating the states, inputs and parameters of a coupled model: 

1. Aerodynamic Model
2. Structural Model
3. Dynamics Model
4. Control Surface Model
5. Controller Model
6. Coupling Model
