# AerostructuralDynamics

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/AerostructuralDynamics.jl/dev)
![](https://github.com/byuflowlab/AerostructuralDynamics.jl/workflows/Run%20tests/badge.svg)

*Aerostructural Dynamics Analysis and Simulation Framework*

Author: Taylor McDonnell

**AerostructuralDynamics** is an extensible multi-fidelity modeling and analysis framework which is designed to simulate the behavior of coupled and/or decoupled aerodynamic, structural, and/or rigid-body dynamics models.  

![](typical-section-flutter-mode.gif)

## Package Features
 - Provides a framework for coupling multiple models together for analysis and/or simulation.
 - Defines a variety of aerodynamic, structural, dynamics, control surface, and control system models
   - Aerodynamics Models:
     - Steady and/or Quasi-Steady Thin Airfoil Theory (2D)
     - Wagner's Function (2D)
     - Peters' Finite State (2D)
     - Lifting Line (3D)
   - Structural Dynamics Models:
     - Two Degree of Freedom Typical Section (2D)
     - Geometrically Exact Beam Theory (As implemented by [GXBeam](https://github.com/byuflowlab/GXBeam.jl)) (3D)
   - Dynamics Models:
     - Rigid Body (3D)
   - Control Surface Models:
     - Steady State Linear Flap (2D)
     - Lifting Line Flaps (3D)
   - Control System Models:
     - Trim (3D)
 - Interfaces with [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
 - Verified and/or validated against theoretical, computational, and/or experimental results (see the [examples](@ref Examples))
 - May be easily extended to include additional models and/or model couplings.  (See the [developer guide](@ref Developer's Guide), pull requests are welcome)

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

## Usage

See the [Getting Started](@ref Getting Started) section of the documentation.