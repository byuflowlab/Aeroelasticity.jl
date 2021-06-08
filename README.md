# AerostructuralDynamics

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/dev)
![](https://github.com/byuflowlab/AerostructuralDynamics.jl/workflows/Run%20tests/badge.svg)

*Aerostructural Dynamics Analysis and Simulation Framework*

Author: Taylor McDonnell

**AerostructuralDynamics** is an analysis tool which simulates the behavior of coupled or decoupled aerodynamic, structural, and/or rigid-body dynamics models.

## Package Features
 - Provides a framework for coupling aerodynamic and structural models for analysis and/or simulation.
 - Interfaces with DifferentialEquations to simulate the behavior of aerostructural models. #TODO
 - Defines a variety of aerodynamic and structural models
   - Aerodynamic Models:
     - Steady (2D)
     - Quasi-Steady (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
     - Steady Vortex Lattice (3D) #TODO
     - Quasi-Steady Vortex Lattice (3D) #TODO
     - Unsteady Vortex Lattice (3D) #TODO
   - Structural Dynamics Models:
     - Typical Section (2D)
     - Geometrically Exact Beam Assemblies (3D) #TODO
   - Dynamics Models:
     - Rigid Body (3D) #TODO
 - Defines the following coupled models
   - Typical Section with:
     - Steady (2D)
     - Quasi-Steady (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
   - Geometrically Exact Beam Assemblies with:
     - Quasi-Steady (2D) #TODO
     - Wagner's Function (2D) #TODO
     - Peters' Finite State (2D) #TODO
     - Steady Vortex Lattice (3D) #TODO
     - Quasi-Steady Vortex Lattice (3D) #TODO
     - Unsteady Vortex Lattice (3D) #TODO
 - Verifications and/or validations of implemented models (see the [examples](https://flow.byu.edu/AerostructuralDynamics.jl/dev))
## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

## Usage

See the [examples](https://flow.byu.edu/AerostructuralDynamics.jl/dev)
