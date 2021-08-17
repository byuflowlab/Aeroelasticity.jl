# Steady Thin Airfoil Theory

![](../airfoil.svg)

## Theory

This model is a steady linear 2D aerodynamic model derived from thin airfoil theory.  As it is a steady model, it has no state variables and/or state equations.  At the reference location, located ``a b`` aft of the semi-chord, the normal force and moment per unit span are defined as
```math
\mathcal{N} = a_0 \rho_\infty u^2 b \alpha_\text{eff} \\
\mathcal{M} = b \left(\frac{1}{2} + a \right) \mathcal{N}
```
where ``u`` is the local freestream velocity in the chordwise direction, ``a`` defines the reference location, ``b`` is the semi-chord, ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density. and ``\alpha_\text{eff}`` is the effective angle of attack.  The effective angle of attack for this model is defined as
```math
\alpha_\text{eff} = -\frac{v}{u} - \alpha_0
```
where ``v`` is the local freestream velocity in the airfoil normal direction and ``\alpha_0`` is the zero lift angle of attack.

## Type Definition

```@docs
Steady
```

## Constructors

```@docs
Steady()
```

## Example Usage
 - [Aeroelastic Analysis of a Typical Section](@ref)
