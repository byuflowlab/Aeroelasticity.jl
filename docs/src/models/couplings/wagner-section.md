# [`Wagner`](@ref) + [`Section`](@ref)

## Type Definition

```@docs
WagnerSection
```

## Example Usage
 - [Aeroelastic Analysis of a Typical Section](@ref section-stability)
 - [Time Domain Simulation of a Typical Section](@ref section-simulation)
 
## Theory

This model is defined by coupling Wagner's function unsteady aerodynamics

![](../../assets/airfoil-drawing.svg)

with the typical section model.  

![](../../assets/section-drawing.svg)

To facilitate this coupling, the freestream velocity components ``u`` and ``v`` are assumed to be aligned with the undeflected chordwise and normal directions, respectively, so that
```math
u \approx U_\infty \\
v \approx \dot{h} \\
\omega \approx \dot{\theta}
```
where ``U_\infty`` is the freestream velocity magnitude, ``\theta`` is pitch, and ``h`` is plunge. To capture the effect of twist on the circulatory lift (since it is no longer implicitly modeled by the ``\frac{v}{u}`` quantity) twist is added to the effective angle of attack from the Wagner's function model so that the effective angle of attack is now given by
```math
\alpha_\text{eff} = \left(\theta - \frac{v}{u} + \frac{b}{u} \left( \frac{1}{2} - a \right) \omega - \alpha_0 \right) \phi(0) + \frac{\lambda_1}{u} + \frac{\lambda_2}{u}
```
The original expression for the effective angle of attack may be used by defining the new variable ``\bar{v} = u \theta + v`` such that
```math
\alpha_\text{eff} = \left(-\frac{\bar{v}}{u} + \frac{b}{u} \left( \frac{1}{2} - a \right) \omega - \alpha_0 \right) \phi(0) + \frac{\lambda_1}{u} + \frac{\lambda_2}{u}
```
A small angle assumption is also used to define the lift about the reference location as
```math
\mathcal{L} \approx \mathcal{N}
```
where ``\mathcal{N}`` is the normal force per unit span at the reference location.

