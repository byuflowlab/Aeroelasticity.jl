# Steady Thin Airfoil Theory

![](../airfoil.svg)

## Type Definition

```@docs
Steady
```

## Constructors

```@docs
Steady()
```

## Theory

This model is two-dimensional aerodynamic model which is based on the results of steady thin airfoil theory.  It does not have any state variables or rate equations.  Aerodynamic forces are assumed to instantaneously reach steady state values when this model is used.

### Airfoil Lift

The lift coefficient of an airfoil according to thin airfoil theory is given by ``a_0 (\alpha - \alpha_0)``, where ``a_0 = 2\pi`` is the airfoil's lift slope, ``\alpha`` is the airfoil's angle of attack, and ``\alpha_0`` is the airfoil's zero lift angle of attack. Dimensionalized, this lift force is
```math
\mathcal{L} = a_0 \rho_\infty b U_\infty^2 (\alpha - \alpha_0)
```
where ``\rho_\infty`` is the freestream air density, ``b`` is the semi-chord, and ``U_\infty`` is the freestream velocity.

### Normal and Axial Force

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

### Pitching Moment

Thin airfoil theory may be used to show that the airfoil quarter chord is the theoretical location of the aerodynamic center.  It may also be used to derive pitching moment coefficients for various airfoil shapes.  To accomodate multiple airfoil shapes, we introduce the quarter-chord moment coefficient ``c_{m_0}`` as an additional parameter.  Using this coefficient, the pitching moment at a location ``a b`` aft of the airfoil mid-chord is defined as
```math
\mathcal{M} = 2 \rho b^2 U_\infty^2 c_{m_0} + b \left(\frac{1}{2} + a \right) \mathcal{N}
```

### Compressibility Correction

At this point, a compressibility correction may be applied to the results of thin airfoil theory in order to extend their applicability.  Applying a Prandtl-Glauert compressibility correction, the normal force, axial force, and pitching moment become
```math
\mathcal{N}_\text{compressible} = \frac{\mathcal{N}}{\sqrt{1 - M^2}} \\
\mathcal{A}_\text{compressible} = \frac{\mathcal{A}}{\sqrt{1 - M^2}} \\
\mathcal{M}_\text{compressible} = \frac{\mathcal{M}}{\sqrt{1 - M^2}}
```
where ``M`` is the Mach number, defined as ``\frac{V_\infty}{c}`` where ``c`` is the air speed of sound. 

### Skin Friction Correction


After the Prandtl-Glauert compressibility correction has been applied, an extra force in the axial direction ``\mathcal{F}_f`` may be added to account for skin friction.  The magnitude of this force is scaled using the ``c_{d_0}`` coefficient.

```math
\mathcal{F}_f = ρ b U_\infty^2 c_{d_0}
```

## Examples
 - [Aeroelastic Analysis of a Typical Section](@ref)
