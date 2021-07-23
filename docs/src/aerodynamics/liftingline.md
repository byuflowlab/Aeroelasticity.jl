# Lifting Line Model

## Theory

Two-dimensional aerodynamic models may be applied in the context of a three-dimensional analysis by applying these models at multiple chordwise sections along the span of one or more lifting surfaces.  This type of model is applicable when spanwise flow effects are negligible, which is often the case for high aspect ratio wings.

The lifting line model implemented in this package assumes that the aerodynamics of each section is independent of the aerodynamics of the other sections, except as coupled through other models.  The state variables and inputs for this model correspond to the state variables, inputs, and parameters of each of the two-dimensional aerodynamic models, concatenated.  Rate equations are also concatenated.  

When coupled with a structural model, aircraft linear and angular accelerations are obtained from the structural model and transformed into the (deformed) local beam frame using an appropriate transformation matrix.  The local freestream velocities/accelerations and pitch rates/accelerations are then defined by a subset of the transformed linear and angular accelerations and cross-flow effects are neglected.   An inverse transformation may then be performed to transform the local aerodynamic forces/moments into the reference frame used by the structural model.

## Type Definition

```@docs
LiftingLine
```

## Constructors

```@docs
LiftingLine(models)
```

## Example Initialization

```@example wagner
using AerostructuralDynamics #hide
model = LiftingLine{4}(Wagner())
nothing #hide
```
