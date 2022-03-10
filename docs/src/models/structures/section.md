# [Typical Section Model](@id typical-section-model)

![](../../assets/section-drawing.svg)

## Type Definition

```@docs
Section
```

## Constructors

```@docs
Section()
```

## Theory

The equations of motion for this model are:
```math
\begin{bmatrix} m & S_\theta \\ S_\theta & I_\theta \end{bmatrix}
\begin{bmatrix} \ddot{h} \\ \ddot{\theta} \end{bmatrix} +
\begin{bmatrix} k_h & 0 \\ 0 & k_h \end{bmatrix}
\begin{bmatrix} h \\ \theta \end{bmatrix} =
\begin{bmatrix} -\mathcal{L} \\ \mathcal{M} \end{bmatrix}
```
where ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, ``I_\theta`` is the mass moment of inertia, ``\mathcal{L}`` is the lift per unit span, and ``\mathcal{M}`` is the moment per unit span.

## Examples
 - [Aeroelastic Analysis of a Typical Section](@ref)
 - [Time Domain Simulation of a Typical Section](@ref)
