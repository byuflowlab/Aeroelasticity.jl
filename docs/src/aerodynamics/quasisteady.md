# Quasi-Steady Thin Airfoil Theory Model

## Theory

This model is a quasi-steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u \dot{\theta} - a b \ddot{\theta} \right) \\
{M'}_{\frac{1}{4}} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right]
```
where ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density, ``u`` is the local freestream velocity in the chordwise direction, ``v`` is the local freestream velocity in the normal direction, ``\theta`` is the pitch angle, ``a`` defines the reference location, ``b`` is the semichord, and ``\alpha_\text{eff}`` is the effective angle of attack.  The reference location is located ``a b`` aft of the semichord.

The effective angle of attack ``\alpha`` for this model is given by
```math
\alpha_\text{eff} = -\frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \dot{\theta} - \alpha_0
```
where ``\alpha_0`` is the zero lift angle of attack.

At the reference location, the lift and moment are
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\dot{\theta} - a b \ddot{\theta} \right) \\
M' = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right] + b \left(\frac{1}{2} + a \right) L'
```

## Documentation

```@docs
QuasiSteady()
```
