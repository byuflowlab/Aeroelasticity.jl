# Steady Thin Airfoil Theory

## Theory

This model is a steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} \\
{M'}_{\frac{1}{4}} = 0
```
where ``a_0`` is the lift curve slope, ``\rho`` is the air density, ``u`` is the local freestream velocity in the chordwise direction, ``b`` is the semichord, and ``\alpha_\text{eff}`` is the effective angle of attack.  For this model, the effective angle of attack is
```math
\alpha_\text{eff} = -\frac{v}{u}
```
where ``v`` is the local freestream velocity in the airfoil normal direction and ``\alpha_0`` is the zero lift angle of attack.

For coupling with structural models, it is convenient to be able to define the lift and moment at an arbitrary reference location.  Defining the reference location to be ``a b`` aft of the semichord, the lift and moment may be expressed as
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} \\
M' = b \left(\frac{1}{2} + a \right) L'
```

## Documentation

```@docs
Steady()
```
