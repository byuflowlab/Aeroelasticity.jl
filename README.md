# Aeroelasticity.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/Aeroelasticity.jl/dev)
![](https://github.com/byuflowlab/Aeroelasticity.jl/workflows/Run%20tests/badge.svg)

*A Modular Multi-Fidelity Aeroelastic Analysis and Simulation Framework*

Author: Taylor McDonnell

**Aeroelasticity.jl** is a modular multi-fidelity aeroelastic analysis and simulation package.  The intent of this package is to facilitate defining and simulating the behavior of both 2D and 3D aeroelastic systems.  Currently, this package only models the aeroelasticity of a two-degree of freedom typical section model, but additional models will be added soon.

![](typical-section-flutter-mode.gif)

## Package Features
 - Facilitates defining and simulating the behavior of both 2D and 3D aeroelastic systems
 - Contains a number of predefined aerodynamic, structural, and aeroelastic models including:
   - Aerodynamic Models 
     - Steady Thin Airfoil Theory (2D)
     - Quasi-Steady Thin Airfoil Theory (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
   - Structural Models
     - Two Degree of Freedom Typical Section Model (2D)
 - Supports multiple analysis types including:
   - Steady State Analyses
   - Eigenvalue Analyses
   - Time-Marching Analyses
 - Directly interfaces with [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
 - Verified/Validated against theoretical, computational, and/or experimental results (see the [examples](@ref Examples))
 - May be easily extended to include additional models and/or model couplings.

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

## Usage

See the [Getting Started](@ref Getting Started) section of the documentation.
