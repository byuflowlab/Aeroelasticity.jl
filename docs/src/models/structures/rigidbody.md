# Rigid Body Model

## Theory

A basic six degree of freedom rigid body model based on the [SixDOF](https://github.com/byuflowlab/SixDOF.jl) package.

## Type Definition

```@docs
RigidBody
```

## Constructors

```@docs
RigidBody()
```

## Example Initialization

```@example rigid-body
using Aeroelasticity #hide
model = RigidBody()
nothing #hide
```
