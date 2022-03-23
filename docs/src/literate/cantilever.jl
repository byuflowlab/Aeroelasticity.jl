# # [Aeroelastic Analysis of a Highly Flexible Cantilever Wing](@id cantilever)
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
# we will use a lifting line model. To model the wing's structure, we will use a 
# geometrically exact beam theory model.

using Aeroelasticity, GXBeam, DifferentialEquations, LinearAlgebra

## discretization
N = 8 # number of elements

## geometric properties
span = 16 # m
chord = 1 # m

## structural section properties
GJ = 1e4 # N*m^2 (torsional rigidity)
EIcc = 2e4 # N*m^2 (flat bending rigidity)
EInn = 4e6 # N*m^2 (chord bending rigidity)
μ = 0.75 # kg/m (mass per unit length)
i11 = 0.1 # kg*m (moment of inertia about elastic axis)
i22 = 0.1*i11 # moment of inertia about beam y-axis
i33 = 0.9*i11 # moment of inertia about beam z-axis

## freestream properties
Vinf = vcat(0.1, 1:35) # m/s (velocity)
α = 2*pi/180 # angle of attack

## aerodynamic section properties
a = 0.0 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
rho = 0.088 # kg/m^3 (air density)
c = 343 # air speed of sound
a0 = 2*pi # lift slope (for each section)
α0 = 0 # zero lift angle of attack (for each section)
cd0 = 0
cm0 = 0

## define geometry
xpt = range(0, 0, length=N+1) # point x-coordinates (in body frame)
ypt = range(0, span, length=N+1) # point y-coordinates (in body frame)
zpt = range(0, 0, length=N+1) # point z-coordinates (in body frame)
points = [[xpt[i],ypt[i],zpt[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
frames = fill([0 1 0; 1 0 0; 0 0 -1], N) # local to body frame transformation
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIcc, 1/EInn]), N) # compliance matrix
mass = fill(Diagonal([μ, μ, μ, i11, i22, i33]), N) # mass matrix
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

prescribed = Dict(
    ## fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

system = System(assembly, false)

## construct section models
section_models = fill(Peters{6}(), N)

## construct aerodynamic model
aerodynamic_model = LiftingLine(section_models)

## construct structural model
structural_model = GXBeamAssembly(assembly, prescribed)

## define coupled model
model = assemble_model(;
    aerodynamic_model = aerodynamic_model,
    structural_model = structural_model)

## current time
t = 0.0

## eigenvalue/eigenvector storage
nev = 12*N
λ = zeros(ComplexF64, nev, length(Vinf))
Uλ = zeros(ComplexF64, nev, number_of_states(model), length(Vinf))
Vλ = zeros(ComplexF64, number_of_states(model), nev, length(Vinf))

## initial guess for state variables
x0 = zeros(number_of_states(model))

## loop through each velocity
for i = 1:length(Vinf)

    ## define parameter vector
    p = assemble_parameters(model;
        aerodynamic_parameters = (; 
            sections = fill((a=a, b=b, a0=a0, alpha0=α0, cd0=cd0, cm0=cm0), N)
            ),
        structural_parameters = (; 
            assembly = assembly
            ),
        additional_parameters = (; 
            v_f = [-Vinf[i]*cos(α), 0, -Vinf[i]*sin(α)],
            rho = rho,
            c = c,
            )
        )

    ## construct ODE function
    f = ODEFunction(model)

    ## find state variables corresponding to steady state operating conditions
    sol = solve(SteadyStateProblem(f, x0, p), SSRootfind())

    ## use state variables from steady state operating conditions
    x = sol.u

    ## linearize about steady state operating conditions
    K, M = linearize(model, x, p)

    ## perform linear stability analysis
    λi, Uλi, Vλi = get_eigen(model, K, M; nev)

    ## correlate eigenvalues
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
    xlim = (0, 35),
    xtick = 0:5:35,
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
    ylim = (-12, 8),
    ytick = -12:4:8,
    ylabel = "Damping Ratio",
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
            real.(λi)./abs.(λi)*100,
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

# # We can visualize the flutter mode using GXBeam's built in interface with WriteVTK

# ## flutter velocity
# Vf = 135

# ## define parameter vector
# p = assemble_parameters(model;
#     aerodynamic_parameters = (; 
#         sections = fill((a=a, b=b, a0=a0, alpha0=α0, cd0=cd0, cm0=cm0), N)
#         ),
#     structural_parameters = (; 
#         assembly = assembly
#         ),
#     additional_parameters = (; 
#         v_f = [-Vinf[i]*cos(α), 0, Vinf[i]*sin(α)],
#         rho = rho,
#         c = c,
#         )
#     )

# ## construct ODE function
# f = ODEFunction(model)

# ## find state variables corresponding to steady state operating conditions
# sol = solve(SteadyStateProblem(f, x0, p), SSRootfind())

# ## use state variables from steady state operating conditions
# x = sol.u

# ## linearize about steady state operating conditions
# K, M = linearize(model, x, p)

# ## perform linear stability analysis
# λf, Uλf, Vλf = get_eigen(model, K, M; nev)

# ## find index corresponding to the flutter mode
# iλ = argmin(abs.(real.(λf)))

# ## find indices corresponding to the GXBeam state variables
# _, igx = state_indices(model)

# ## flutter eigenvalue
# eigenvalue = 0 + imag(λf[iλ])*1im

# ## post-process GXBeam state variables
# state = AssemblyState(system, assembly, xf[igx]; prescribed_conditions=prescribed)

# ## post-process GXBeam eigenvector variables
# eigenstate = AssemblyState(system, assembly, Vλf[igx,iλ]; prescribed_conditions=prescribed)

# ## create a 4% thick parabolic airfoil
# f = (x) -> 0.16*x*(1 - x)
# x = 0:0.2:1.0
# zu = f.(x)
# zl = -zu
# x = chord .* vcat(x, x[end-1:-1:1])
# y = chord .* zeros(length(x))
# z = chord .* vcat(zu, zl[end-1:-1:1])
# sections = vcat(x', y', z')

# ## write the response to vtk files for visualization using ParaView
# write_vtk("cantilever-flutter-mode", assembly, state, eigenvalue, eigenstate; sections=sections,
#     mode_scaling=100)

# #md # ![](../assets/goland-flutter-mode.gif)
