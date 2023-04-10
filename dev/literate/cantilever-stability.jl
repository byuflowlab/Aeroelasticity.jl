# # [Aeroelastic Analysis of a Highly Flexible Wing](@id cantilever)
#
# In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis
# of a highly flexible cantilever wing.
#
# ![](../assets/cantilever-wing.png)
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`cantilever.ipynb`](@__NBVIEWER_ROOT_URL__/examples/cantilever.ipynb).
#-
#
# The wing we are considering in this example was created by modifying Daedalus aircraft
# data and is therefore representative of a high-altitude long-endurance wing. It has a
# 16 meter span (from root to tip) and a 1 meter chord. To model the wing's aerodynamics,
# we use a lifting line model. To model the wing's structure, we use a geometrically exact
# beam theory model.

using Aeroelasticity, GXBeam, DifferentialEquations, LinearAlgebra

## --- Initial Setup --- #

## discretization
N = 8 # number of elements

## geometric properties
span = 16 # m
chord = 1 # m
xref = 0.5 # normalized reference location (from leading edge)
xcg = 0.5 # center of gravity (from leading edge)

## stiffness properties
GJ = 1e4 # N*m^2 (torsional rigidity)
EIyy = 2e4 # N*m^2 (flat bending rigidity)
EIzz = 4e6 # N*m^2 (chord bending rigidity)

# inertial properties
mu = 0.75 # kg/m (mass per unit length)
i11 = 0.1 # kg*m (moment of inertia about elastic axis)
i22 = 0.0375 # moment of inertia about beam y-axis
i33 = 0.0625 # moment of inertia about beam z-axis

## freestream properties
Vinf = 10.0 # m/s (velocity)
rho = 0.0889 # kg/m^3 (air density at 20 km)
c = 343 # m/s (air speed of sound)
alpha = 2*pi/180 # angle of attack

## aerodynamic section properties
a = xref - 0.5 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
a0 = 2*pi # lift slope (for each section)
alpha0 = 0 # zero lift angle of attack (for each section)
cd0 = 0
cm0 = 0

## define geometry (assume NED coordinate frame)
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
damping = fill(fill(1e-4, 6), N) # stiffness proportional structural damping coefficients
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass, damping)

prescribed_conditions = Dict(
    ## fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
)

## define GXBeam system
system = DynamicSystem(assembly)

## --- Define Submodels --- #

## construct section models
section_models = fill(Peters{6}(), N)

## construct aerodynamic model
aerodynamic_model = LiftingLine(section_models)

## construct structural model
structural_model = GXBeamAssembly(system; structural_damping=true)

## define submodels
submodels = (aerodynamic_model, structural_model)

# --- Define Initial Parameters --- #

V = [-Vinf*cos(alpha), 0.0, -Vinf*sin(alpha)] # m/s (freestream velocity)

## define parameters for each lifting line section
section_parameters = fill([a, b, a0, alpha0, cd0, cm0], N)

## define parameters for the lifting line model
liftingline_parameters = LiftingLineParameters(section_parameters)

## define parameters for the geometrically exact beam theory model
gxbeam_parameters = GXBeamParameters(assembly)

# define parameters for the coupling
coupling_parameters = LiftingLineGXBeamParameters(V, rho, c;
    prescribed_conditions = prescribed_conditions,
    gravity = [-9.81*sin(alpha), 0, 9.81*cos(alpha)])

## combine parameters
parameters = (liftingline_parameters, gxbeam_parameters, coupling_parameters)

## --- Define Coupled Model --- #

model = CoupledModel(submodels, parameters; symbolic=false)

## --- Perform Analysis --- #

## loop through freestream velocities
Vinf = vcat(0.1, 0.25:0.25:35)

## eigenvalue/eigenvector storage
nev = 12*N
λ = zeros(ComplexF64, nev, length(Vinf))
Uλ = zeros(ComplexF64, nev, number_of_states(model), length(Vinf))
Vλ = zeros(ComplexF64, number_of_states(model), nev, length(Vinf))

## initial guess for state variables
x0 = zeros(number_of_states(model))

## loop through each velocity
for i = 1:length(Vinf)

    ## --- Update Parameters --- #

    V = [-Vinf[i]*cos(alpha), 0.0, -Vinf[i]*sin(alpha)] # m/s (freestream velocity)

    ## define parameters for each lifting line section
    section_parameters = fill([a, b, a0, alpha0, cd0, cm0], N)

    ## define parameters for the lifting line model
    liftingline_parameters = LiftingLineParameters(section_parameters)

    ## define parameters for the geometrically exact beam theory model
    gxbeam_parameters = GXBeamParameters(assembly)

    ## define parameters for the coupling
    coupling_parameters = LiftingLineGXBeamParameters(V, rho, c;
        prescribed_conditions = prescribed_conditions,
        gravity = [-9.81*sin(alpha), 0, 9.81*cos(alpha)])

    ## combine parameters
    parameters = (liftingline_parameters, gxbeam_parameters, coupling_parameters)

    ## --- Perform Analysis --- #

    ## define an ODEFunction for this model
    f = ODEFunction(model, parameters)

    ## find equilibrium point
    sol = solve(NonlinearProblem(SteadyStateProblem(f, x0, parameters)))

    ## use state variables from steady state operating conditions
    x = sol.u

    ## linearize about steady state operating conditions
    K, M = linearize(model, x, parameters)

    ## perform linear stability analysis
    λi, Uλi, Vλi = sparse_eigen(K, M; nev=nev)

    ## --- Correlate Eigenvalues --- #

    if i > 1
        ## previous left eigenvector matrix
        Uλpi = Uλ[:,:,i-1]

        ## use correlation matrix to correlate eigenmodes
        perm, corruption = Aeroelasticity.correlate_eigenmodes(Uλpi, M, Vλi)

        ## re-arrange eigenmodes
        λi = λi[perm]
        Uλi = Uλi[perm,:]
        Vλi = Vλi[:,perm]
    end

    ## save eigenvalues/eigenvectors
    λ[:,i] = λi
    Uλ[:,:,i] = Uλi
    Vλ[:,:,i] = Vλi

    ## update initial guess for the state variables
    x0 .= x
end

# To identify the flutter speed and frequency, we can plot the results.

using Plots
pyplot()

sp1 = plot(
    xlim = (0, 50),
    xtick = 0:5:50,
    xlabel = "Velocity (m/s)",
    ylim = (0, 350),
    ytick = 0:50:350,
    ylabel = "Frequency (rad/s)",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topright,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true)

sp2 = plot(
    xlim = (0, 35),
    xtick = 0:5:35,
    xlabel = "Velocity (m/s)",
    ylim = (-4, 2),
    ytick = -4:2:2,
    ylabel = "Damping (1/s)",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true)

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(abs.(λi) .<= 1e4)
        plot!(sp1, Vi, imag.(λi),
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end

end

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(abs.(λi) .<= 1e4)
        plot!(sp2, Vi,
            real.(λi),#./sqrt.(real.(λi).^2 + abs.(λi).^2)*100,
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

#jl plot!(show=true)
#md savefig("../assets/cantilever-stability.svg")
#md nothing

#md # ![]("../assets/cantilever-stability.svg")

#-

# These results are similar to those presented by Patil, Hodges, and Cesnik in "Nonlinear
# Aeroelasticity and Flight Dynamics of High Altitude Long-Endurance Aircraft" and
# therefore serve as a verification case for this coupled model.
