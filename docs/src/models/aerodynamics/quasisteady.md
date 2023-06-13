# [Quasi-Steady Thin Airfoil Theory](@id quasi-steady-theory)

![](../../assets/airfoil-drawing.svg)

## Type Definition

```@docs
QuasiSteady
QuasiSteady()
```

## Examples
 - [Aeroelastic Analysis of a Typical Section](@ref section-stability)
 - [Time Domain Simulation of a Typical Section](@ref section-simulation)
 - [Aeroelastic Analysis of the Goland Wing](@ref goland-stability)
 - [Steady State Aeroelastic Analysis of a Highly Flexible Wing](@ref cantilever-steady)
 - [Aeroelastic Stability Analysis of a Highly Flexible Wing](@ref cantilever-stability)

## Theory

This model is a two-dimensional quasi-steady aerodynamic model derived from thin airfoil theory.  It does not have any state variables or rate equations.  Since this model includes no state variables to model wake development, this model assumes that an airfoil's wake reaches steady-state operating conditions instantaneously in response to changes in freestream velocities.

### Normal Force, Axial Force, and Pitching Moment

The combined forces and moments ``a b`` aft of the mid-chord for this model are
```math
\mathcal{N} = a_0 \rho_\infty b u^2 \alpha_\text{eff} + \pi \rho b^2 \left(\dot{v} + u\omega - a b \dot{\omega} \right) \\
\mathcal{A} = -a_0 \rho_\infty b u v \alpha_\text{eff} \\
\mathcal{M} = 2 \rho b^2 u^2 c_{m_0} -\pi \rho_\infty b^3 \left[\frac{1}{2}\dot{v} + u\omega + b \left( \frac{1}{8} - \frac{a}{2} \right) \dot{\omega} \right] + b \left(\frac{1}{2} + a \right) \mathcal{N}
```
where ``u`` is the local freestream velocity in the chordwise direction, ``v`` is the local freestream velocity in the normal direction, ``\omega`` is the freestream angular velocity, ``a`` defines the reference location, ``b`` is the semichord, ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density, and ``\alpha_\text{eff}`` is the effective angle of attack.  The effective angle of attack for this model is defined as
```math
\alpha_\text{eff} = \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega - \alpha_0
```
where ``\alpha_0`` is the zero lift angle of attack.

### Compressibility Correction

At this point, a compressibility correction may be applied to the results of thin airfoil theory in order to extend their applicability.  Applying a Prandtl-Glauert compressibility correction, the normal force, axial force, and pitching moment become
```math
\mathcal{N}_\text{compressible} = \frac{\mathcal{N}}{\beta} \\
\mathcal{A}_\text{compressible} = \frac{\mathcal{A}}{\beta} \\
\mathcal{M}_\text{compressible} = \frac{\mathcal{M}}{\beta}
```
where ``\beta=\sqrt{1 - M^2}`` is the Prandtl-Glauert factor, which is a function of the local section Mach number ``M``.

### Viscous Forces

After the Prandtl-Glauert compressibility correction has been applied, an extra force in the axial direction ``\mathcal{F}_v`` may be added to account for viscous forces.  The magnitude of this force is scaled using the ``c_{d_0}`` coefficient.

```math
\mathcal{F}_v = ρ b u^2 c_{d_0}
```

