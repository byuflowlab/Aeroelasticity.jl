# [`Steady`](@ref) + [`Section`](@ref)

## Type Definition

```@docs
SteadySection
```

## Example Usage
 - [Aeroelastic Analysis of a Typical Section](@ref section-stability)
 - [Time Domain Simulation of a Typical Section](@ref section-simulation) 

## Theory

This model is defined by coupling steady thin airfoil theory aerodynamics

![](../../assets/airfoil-drawing.svg)

with the typical section model.  

![](../../assets/section-drawing.svg)

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


