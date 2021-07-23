# [`Steady`](@ref) + [`TypicalSection`](@ref)

![](../airfoil.svg)
![](../typical-section.svg)

## Theory

This model is defined by coupling steady thin airfoil theory aerodynamics with the typical section model.  By restricting the analysis to small angles, this model is able to assume
```math
u \approx U_\infty \\
v \approx U_\infty \theta \\
\mathcal{L} \approx \mathcal{N}
```
where ``u`` is the chordwise freestream velocity, ``v`` is the normal freestream velocity, ``U_\infty`` is the freestream velocity magnitude, ``\theta`` is the pitch angle, ``\mathcal{L}`` is the lift per unit span and ``\mathcal{N}`` is the normal force per unit span.

## Constructors

```@docs
couple_models(aero::Steady, stru::TypicalSection)
```

## Example Initialization

```@example steady-section
using AerostructuralDynamics #hide
model = couple_models(Steady(), TypicalSection())
nothing #hide
```
