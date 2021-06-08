# Models

## Aerodynamic Models

### Two-Dimensional Models

#### Steady [`Steady`](@ref)

#### Quasi-Steady [`QuasiSteady`](@ref)

#### Wagner's Function [`Wagner`](@ref)

#### Peter's Finite State Model [`Peters`](@ref)

### Three-Dimensional Models

#### Lifting Line [`LiftingLine`](@ref)

#### Vortex Lattice Method [`VLM`](@ref)

## Structural Models

### Two-Dimensional Models

#### Typical Section [`TypicalSection`](@ref)

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M_{\frac{1}{4}} + b \left( \frac{1}{2} + a \right) L
```
where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M_\frac{1}{4}`` is the quarter-chord moment per unit span.

### Three-Dimensional Models

#### Rigid Body [`RigidBody`](@ref)

#### Geometrically Exact Beam Theory [`GEBT`](@ref)
