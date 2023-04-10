using Aeroelasticity, GXBeam, NLsolve, LineSearches, LinearAlgebra

# --- Initial Setup --- #

# discretization
N = 8 # number of elements

# geometric properties
span = 16 # m
chord = 1 # m
xref = 0.5 # normalized reference location (from leading edge)
xcg = 0.5 # center of gravity (from leading edge)

# stiffness properties
GJ = 1e4 # N*m^2 (torsional rigidity)
EIyy = 2e4 # N*m^2 (flat bending rigidity)
EIzz = 4e6 # N*m^2 (chord bending rigidity)

mu = 0.75 # kg/m (mass per unit length)
i11 = 0.1 # kg*m (moment of inertia about elastic axis)
i22 = 0.0375 # moment of inertia about beam y-axis
i33 = 0.0625 # moment of inertia about beam z-axis

# freestream properties
Vinf = 25.0 # m/s (velocity)
rho = 0.0889 # kg/m^3 (air density)
c = 343 # m/s (air speed of sound)
alpha = 10*pi/180 # angle of attack

# aerodynamic section properties
a = xref - 0.5 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
a0 = 2*pi # lift slope (for each section)
alpha0 = 0.0 # zero lift angle of attack (for each section)
cd0 = 0.01
cm0 = 0

# define geometry (assume NED coordinate frame)
xpt = range(0, 0, length=N+1) # point x-coordinates (in body frame)
ypt = range(0, span, length=N+1) # point y-coordinates (in body frame)
zpt = range(0, 0, length=N+1) # point z-coordinates (in body frame)
points = [[xpt[i],ypt[i],zpt[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
e1 = [0, 1,  0] # beam x-axis
e2 = [1, 0,  0] # beam y-axis
e3 = [0, 0, -1] # beam z-axis
frames = fill([e1 e2 e3], N) # local to body frame transformation
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIyy, 1/EIzz]), N) # compliance matrix
mass = fill(Diagonal([mu, mu, mu, i11, i22, i33]), N) # mass matrix
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

prescribed_conditions = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
)

# define GXBeam system
system = DynamicSystem(assembly)

# --- Define Submodels --- #

# construct section models
section_models = fill(Peters{6}(), N)

# construct aerodynamic model
aerodynamic_model = LiftingLine(section_models)

# construct structural model
structural_model = GXBeamAssembly(system; structural_damping=false)

# define submodels
submodels = (aerodynamic_model, structural_model)

# --- Define Initial Parameters --- #

V = [-Vinf*cos(alpha), 0.0, -Vinf*sin(alpha)] # m/s (freestream velocity)

# define parameters for each lifting line section
section_parameters = fill([a, b, a0, alpha0, cd0, cm0], N)

# define parameters for the lifting line model
liftingline_parameters = LiftingLineParameters(section_parameters)

# define parameters for the geometrically exact beam theory model
gxbeam_parameters = GXBeamParameters(assembly)

coupling_parameters = LiftingLineGXBeamParameters(V, rho, c;
    prescribed_conditions = prescribed_conditions,
    gravity = [0, 0, 9.81])

# combine parameters
parameters = (liftingline_parameters, gxbeam_parameters, coupling_parameters)

# --- Define Coupled Model --- #

model = CoupledModel(submodels, parameters; symbolic=false)

# --- Perform Analysis --- #

# state rates equal to zero
dx = zeros(number_of_states(model))

# use previously defined parameters
p = parameters

# set time to zero
t = 0.0

# define residual function
f! = (resid, x) -> residual!(resid, model, dx, x, p, t)

# define jacobian function
j! = (jacob, x) -> state_jacobian!(jacob, model, dx, x, p, t)

# define initial guess
x0 = zeros(number_of_states(model))

# find steady state operating conditions
result = nlsolve(f!, j!, x0; method=:newton, linesearch=LineSearches.BackTracking())

# separate state variables into aerodynamic and structural states
λ, x = separate_states(result.zero, model)

state = AssemblyState(x, system, assembly; prescribed_conditions)

write_vtk("cantilever-steady", assembly, state)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

