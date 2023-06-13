# [Steady Thin Airfoil Theory](@id steady-theory)

![](../../assets/airfoil-drawing.svg)

## Type Definition

```@docs
Steady
Steady()
```

## Examples
 - [Aeroelastic Analysis of a Typical Section](@ref section-stability)
 - [Time Domain Simulation of a Typical Section](@ref section-simulation)
 - [Aeroelastic Analysis of the Goland Wing](@ref goland-stability)
 - [Steady State Aeroelastic Analysis of a Highly Flexible Wing](@ref cantilever-steady)
 - [Aeroelastic Stability Analysis of a Highly Flexible Wing](@ref cantilever-stability)

## Theory

This model is a two-dimensional aerodynamic model which is based on the results of steady thin airfoil theory.  As it is a steady state model, it does not have any state varables or rate equations.  In other words, this model assumes that aerodynamic forces instantaneously assume their steady state values.

## Normal and Axial Forces

The lift coefficient of an airfoil according to thin airfoil theory is given by ``a_0 (\alpha - \alpha_0)``, where ``a_0 = 2\pi`` is the airfoil's lift slope, ``\alpha`` is the airfoil's angle of attack, and ``\alpha_0`` is the airfoil's zero lift angle of attack. Dimensionalized, this lift force is
```math
\mathcal{L} = a_0 \rho_\infty b U_\infty^2 (\alpha - \alpha_0)
```
where ``\rho_\infty`` is the freestream air density, ``b`` is the semi-chord, and ``U_\infty`` is the freestream velocity.

The lift force may be expressed in terms of normal and axial components as
```math
\mathcal{N} = a_0 \rho_\infty b U_\infty^2 (\alpha - \alpha_0) cos(\alpha) \\
\mathcal{A} = -a_0 \rho_\infty b U_\infty^2 (\alpha - \alpha_0) sin(\alpha)
```
Using the substitutions
```math
cos(\alpha) = \frac{u}{U_\infty} \quad sin(\alpha) = \frac{v}{U_\infty}
```
where ``u`` is the tangential velocity and ``v`` is the normal velocity leads to the following expressions for the normal and axial forces
```math
\mathcal{N} = a_0 \rho_\infty b U_\infty u (\alpha - \alpha_0) \\
\mathcal{A} = -a_0 \rho_\infty b U_\infty v (\alpha - \alpha_0)
```
Using a small angle approximation allows us to assume ``\alpha \approx sin(\alpha) = \frac{v}{U_\infty} ``, which allows us to reduce our expressions for the normal and axial force to
```math
\mathcal{N} = a_0 \rho_\infty b u (v - U_\infty \alpha_0) \\
\mathcal{A} = -a_0 \rho_\infty b v (v - U_\infty \alpha_0)
```
If we further assume that ``u >> v`` (which is reasonable considering our small angle assumption) then
```math
U_\infty = \sqrt{u^2 + v^2} \approx = u
```
and the expressions for the normal and axial force further reduce to
```math
\mathcal{N} = a_0 \rho_\infty b u^2 \alpha_\text{eff} \\
\mathcal{A} = -a_0 \rho_\infty b u v \alpha_\text{eff}
```
where ``\alpha_\text{eff}`` is the effective angle of attack, defined for this model as
```math
\alpha_\text{eff} = \frac{v}{u} - \alpha_0
```

## Pitching Moment

Thin airfoil theory may be used to show that the airfoil quarter chord is the theoretical location of the aerodynamic center.  It may also be used to derive pitching moment coefficients for various airfoil shapes.  To accomodate multiple airfoil shapes, we introduce the quarter-chord moment coefficient ``c_{m_0}`` as an additional parameter.  Using this coefficient, the pitching moment at a location ``a b`` aft of the airfoil mid-chord may be defined as
```math
\mathcal{M} = 2 \rho b^2 u^2 c_{m_0} + b \left(\frac{1}{2} + a \right) \mathcal{N}
```

## Compressibility Correction

At this point, a compressibility correction may be applied to the results of thin airfoil theory in order to extend their applicability.  Applying a Prandtl-Glauert compressibility correction, the normal force, axial force, and pitching moment become
```math
\mathcal{N}_\text{compressible} = \frac{\mathcal{N}}{\beta} \\
\mathcal{A}_\text{compressible} = \frac{\mathcal{A}}{\beta} \\
\mathcal{M}_\text{compressible} = \frac{\mathcal{M}}{\beta}
```
where ``\beta=\sqrt{1 - M^2}`` is the Prandtl-Glauert factor, which is a function of the local section Mach number ``M``.

## Viscous Forces

After the Prandtl-Glauert compressibility correction has been applied, an extra force in the axial direction ``\mathcal{F}_v`` may be added to account for viscous forces.  The magnitude of this force is scaled using the ``c_{d_0}`` coefficient.

```math
\mathcal{F}_v = ρ b u^2 c_{d_0}
```


