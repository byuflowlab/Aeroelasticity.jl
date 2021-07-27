# [`Steady`](@ref) + [`TypicalSection`](@ref)

## Theory

This model is defined by coupling steady thin airfoil theory aerodynamics

![](../airfoil.svg)

with the typical section model.  

![](../typical-section.svg)

By making use of a small angle assumption, the freestream velocity components are defined as
```math
u \approx U_\infty \\
v \approx U_\infty \theta \\
\omega \approx 0
```
where ``U_\infty`` is the freestream velocity magnitude, ``\theta`` is the pitch angle, ``u`` is the chordwise freestream velocity, ``v`` is the normal freestream velocity, and ``\omega`` is the freestream angular velocity. A small angle assumption is also used to define the lift about the reference location as
```math
\mathcal{L} \approx \mathcal{N}
```
where ``\mathcal{N}`` is the normal force per unit span at the reference location.

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
