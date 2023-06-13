# Aeroelasticity.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/Aeroelasticity.jl/dev)
![](https://github.com/byuflowlab/Aeroelasticity.jl/workflows/Run%20tests/badge.svg)

*A Modular Multi-Fidelity Aeroelastic Analysis and Simulation Framework*

Author: Taylor McDonnell

**Aeroelasticity.jl** is a modular multi-fidelity aeroelastic analysis and simulation package.  The intent of this package is to facilitate defining and simulating the behavior of both 2D and 3D aeroelastic systems, though it may also be used to model other coupled systems.

![](assets/section-simulation.gif)

## Package Features
 - Facilitates defining and simulating the behavior of both 2D and 3D aeroelastic systems
 - Contains a number of predefined aerodynamic, structural, and aeroelastic models including:
   - Aerodynamic Models
     - Steady Thin Airfoil Theory (2D)
     - Quasi-Steady Thin Airfoil Theory (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
     - Lifting Line (3D)
   - Structural Models
     - Two Degree of Freedom Typical Section Model (2D)
     - Rigid Body Model (3D)
     - Geometrically Exact Beam Theory Model (3D)
 - Supports multiple analysis types including:
   - Steady State Analyses
   - Eigenvalue Analyses
   - Time-Marching Analyses (using [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl))
 - Several verification cases for built-in models (see the examples)
 - Provides a modular framework for constructing general monolithic coupled systems

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://github.com/byuflowlab/Aeroelasticity.jl
```

## Usage

See the [Getting Started](@ref guide) section of the documentation.
