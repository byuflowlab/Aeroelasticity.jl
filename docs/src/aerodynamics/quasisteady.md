# [Quasi-Steady Thin Airfoil Theory Model](@id quasi-steady-thin-airfoil-theory)

![](../airfoil.svg)

## Theory

This model is a quasi-steady linear 2D aerodynamic model derived from thin airfoil theory.  As it is a quasi-steady model, it has no state variables and/or state equations.  At the reference location, located ``a b`` aft of the semi-chord, the normal force and moment per unit span are defined as
```math
\mathcal{N} = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\omega - a b \dot{\omega} \right) \\
\mathcal{M} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\omega + b \left( \frac{1}{8} - \frac{a}{2} \right) \dot{\omega} \right] + b \left(\frac{1}{2} + a \right) \mathcal{N}
```
where ``u`` is the local freestream velocity in the chordwise direction, ``v`` is the local freestream velocity in the normal direction, ``\omega`` is the freestream angular velocity, ``a`` defines the reference location, ``b`` is the semichord, ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density, and ``\alpha_\text{eff}`` is the effective angle of attack.  The effective angle of attack for this model is defined as
```math
\alpha_\text{eff} = - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega - \alpha_0
```
where ``\alpha_0`` is the zero lift angle of attack.

## Type Definition

```@docs
QuasiSteady
```

## Constructors

```@docs
QuasiSteady()
```

## Example Usage
 - [Aeroelastic Analysis of a Typical Section](@ref)
