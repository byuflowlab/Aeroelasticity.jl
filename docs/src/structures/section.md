# [Typical Section Model](@id typical-section-model)

![](typical-section.svg)

## Theory

The equations of motion for this model are:
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -\mathcal{L} \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = \mathcal{M}
```
where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M`` is the moment per unit span about the reference point.

## Documentation

```@docs
TypicalSection()
```
