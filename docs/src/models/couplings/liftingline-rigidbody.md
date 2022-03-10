# [`LiftingLine`](@ref) + [`RigidBody`](@ref)

## Theory

This model is defined by coupling the lifting line aerodynamics model with the rigid body model.

## Constructors

```@docs
couple_models(aero::LiftingLine, stru::RigidBody)
```

## Example Initialization

```@example liftingline-rigidbody
using AerostructuralDynamics #hide
model = couple_models(LiftingLine{6}(Wagner()), RigidBody())
nothing #hide
```
