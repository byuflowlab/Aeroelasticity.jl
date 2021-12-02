# Geometrically Exact Beam Theory Model

## Theory

A geometrically exact beam model implemented using the [GXBeam](https://github.com/byuflowlab/GXBeam.jl) package.

## Type Definition

"""
GEBT

Construct a geometrically exact beam theory model, as implemented by the GXBeam package.
State variables are as defined by GXBeam.  Inputs correspond to the external
forces ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or
displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i},
\\theta_{z,i}`` applied to each node, followed by the distributed loads
``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}`` applied to each beam
element, followed by the properties of point masses attached to each beam element ``m, p, 
I_{11}, I_{22}, I_{33}, I_{12}, I_{13}, I_{23}``, followed by the linear and angular 
velocity and accelearation of the system. Parameters correspond to the location ``p_{x}, 
p_{y}, p_{z}`` of each node followed by each beam element's properties. Each beam element's 
properties are defined by a triad which defines the orientation of the beam element 
``e_{1,x}, e_{1,y}, e_{1,z}, e_{2,x}, e_{2,y}, e_{2,z}, e_{3,x}, e_{3,y}, e_{3,z}``, 
followed by the 21 independent entries of the compliance matrix ``C_{11}, C_{12}, C_{13}, 
C_{14}, C_{15}, C_{16}, C_{22}, C_{23}, C_{24}, C_{25}, C_{26}, C_{33}, C_{34}, C_{35}, 
C_{36}, C_{44}, C_{45}, C_{46}, C_{55}, C_{56}, C_{66}``, followed by the beam element's 
inertial properties ``\\mu, x_{m,2}, x_{m,3}, i_{22}, i_{33}, i_{23}``.
"""
GEBT

```@docs
GEBT
```

## Constructors

```@docs
GEBT(assembly, prescribed)
```

## Example Initialization

```@example gxbeam
using AerostructuralDynamics, GXBeam, LinearAlgebra

# discretization
N = 8 # number of elements

# geometric properties
span = 6.096 # m (wing half span)
chord = 1.8288 # m (chord)

# structural section properties
xea = 0.33*chord # m (elastic axis, from leading edge)
EIcc = 9.77e6 # N*m^2 (flat bending rigidity)
GJ = 0.99e6 # N*m^2 (torsional rigidity)
μ = 35.71 # kg/m (mass per unit length)
xcm = 0.43*chord # m (center of mass, from leading edge)
i11 = 8.64 # kg*m (moment of inertia about elastic axis)
i22 = 0.1*i11 # moment of inertia about beam y-axis
i33 = 0.9*i11 # moment of inertia about beam z-axis

# define geometry
xpt = range(0, 0, length=N+1) # point x-coordinates (in body frame)
ypt = range(0, span, length=N+1) # point y-coordinates (in body frame)
zpt = range(0, 0, length=N+1) # point z-coordinates (in body frame)
points = [[xpt[i],ypt[i],zpt[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
frames = fill([0 1 0; 1 0 0; 0 0 -1], N) # local to body frame transformation
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIcc, 0]), N) # compliance matrix
xm2 = xea - xcm
mass = fill([
    μ 0 0 0 0 -μ*xm2;
    0 μ 0 0 0 0;
    0 0 μ μ*xm2 0 0;
    0 0 μ*xm2 i11 0 0;
    0 0 0 0 i22 0;
    -μ*xm2 0 0 0 0 i33], N) # mass matrix
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

# boundary condition initialization
prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

model = GEBT(assembly, prescribed)

nothing #hide
```
