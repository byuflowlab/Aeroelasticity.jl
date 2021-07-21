# AerostructuralDynamics

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/dev)
![](https://github.com/byuflowlab/AerostructuralDynamics.jl/workflows/Run%20tests/badge.svg)

*Aerostructural Dynamics Analysis and Simulation Framework*

Author: Taylor McDonnell

**AerostructuralDynamics** is a multi-fidelity modeling and analysis framework which is designed to simulate the behavior of coupled and/or decoupled aerodynamic, structural, and/or rigid-body dynamics models.

## Package Features
 - Provides a framework for coupling multiple models together for analysis and/or simulation.
 - Defines a variety of aerodynamic, structural, and dynamics models
   - Aerodynamic Models:
     - Steady and/or Quasi-Steady Thin Airfoil Theory (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
     - Lifting Line (3D)
   - Structural Dynamics Models:
     - Two Degree of Freedom Typical Section (2D)
     - Geometrically Exact Beam Theory (As implemented by [GXBeam](https://github.com/byuflowlab/GXBeam.jl)) (3D)
   - Dynamics Models:
     - Rigid Body (3D)
 - Verified and/or validated against theoretical, computational, and/or experimental results (see the [examples](@ref Examples))

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

## Usage

See the [examples](@ref Examples)
