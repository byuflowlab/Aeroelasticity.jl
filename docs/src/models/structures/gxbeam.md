# Geometrically Exact Beam Element Assembly

## Type Definition

```@docs
GXBeamAssembly
GXBeamAssembly()
GXBeamInputs
GXBeamParameters
```

## Examples
 - [Aeroelastic Analysis of the Goland Wing](@ref goland-stability)
 - [Steady State Aeroelastic Analysis of a Highly Flexible Wing](@ref cantilever-steady)
 - [Aeroelastic Stability Analysis of a Highly Flexible Wing](@ref cantilever-stability)
 
## Theory

This model uses geometrically exact beam theory (as implemented by the [GXBeam](https://github.com/byuflowlab/GXBeam.jl) package to model an interconnected assembly of nonlinear beams.  For more details, refer to the documentation for the [GXBeam](https://github.com/byuflowlab/GXBeam.jl) package.

