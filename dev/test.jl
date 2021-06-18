using AerostructuralDynamics, GXBeam, NLsolve, LinearAlgebra

# discretization
N = 8 # number of elements

# geometric properties
span = 16 # m
chord = 1 # m (chord)
xref = 0.5 # normalized reference location (from leading edge)
xcg = 0.5 # center of gravity (from leading edge)

# freestream properties
Vinf = range(1, 30, length=50) # m/s (velocity)
α = 2*pi/180 # angle of attack

# aerodynamic section properties
a = xref - 0.5 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
ρ = 0.088 # kg/m^3 (air density)
a0 = 2*pi # lift slope (for each section)
α0 = 0 # zero lift angle of attack (for each section)

# structural section properties
EIcc = 2e4 # N*m^2 (flat bending rigidity)
EInn = 4e6 # N*m^2 (chord bending rigidity)
GJ = 1e4 # N*m^2 (torsional rigidity)
μ = 0.75 # kg/m (mass per unit span)
i11 = 0.1 # kg*m (rotational inertia per unit span)
i22 = 0.0375 # moment of inertia about beam y-axis
i33 = 0.0625 # moment of inertia about beam z-axis

# define geometry
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

# boundary condition initialization
prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

# distributed load initialization
distributed = Dict()
for i = 1:N
    # distributed load on each beam element
    distributed[i] = GXBeam.DistributedLoads(assembly, i)
end

# structural system initialization
system = GXBeam.System(assembly, keys(prescribed), false)

# construct aerodynamic model
aerodynamic_model = LiftingLine{N}(Peters{6}())

# construct structural model
structural_model = GEBT(system, assembly, prescribed, distributed)

# define simulation models
models = (aerodynamic_model, structural_model)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(models), length(Vinf))

# state variable initial guess
u0 = zeros(number_of_states(models))

# loop through each velocity
for i = 1:length(Vinf)

    println("Vinf: ", Vinf[i])

    # set parameters and current time
    p_aero = vcat(fill([a, b, a0, α0], N)...)
    p_stru = Float64[]
    p_additional = [-Vinf[i]*cos(α), 0, -Vinf[i]*sin(α), ρ]
    p = vcat(p_aero, p_stru, p_additional)
    t = 0

    # find state variables corresponding to steady state operating conditions
    fresid = u -> get_rates(models, u, get_inputs(models, u, p, t), p, t)
    sol = nlsolve(fresid, u0)
    u = sol.zero

    # calculate the inputs corresponding to steady state operating conditions
    y = get_inputs(models, u, p, t)

    # calculate the mass matrix corresponding to steady state operating conditions
    M = get_mass_matrix(models, u, y, p, t)

    # calculate the jacobian corresponding to steady state operating conditions
    J = get_state_jacobian(models, u, y, p, t)

    # solve the generalized eigenvalue problem
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

    # update initial guess for the state variables
    u0 .= u
end

using Plots
pyplot()

default(
    titlefontsize = 14,
    legendfontsize = 11,
    guidefontsize = 14,
    tickfontsize = 11,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true,
    framestyle = :zerolines)

sp1 = plot(
    xlim = (1, 30),
    xtick = vcat(1, 5:5:30),
    xlabel = "Velocity (m/s)",
    ylim = (0, 90),
    ytick = 0.0:10:90,
    ylabel = "Frequency (rad/s)",
    legend = :topright
    )

sp2 = plot(
    xlim = (1, 30),
    xtick = vcat(1, 5:5:30),
    xlabel = "Velocity (m/s)",
    ylim = (-12, 8),
    ytick = -12:4:8,
    ylabel = "Damping Ratio %",
    legend = :topleft
    )

for i = 1:size(λ, 1)

    idx = findall(x -> abs(x) < 500, λ[i,:])
    Vi = Vinf[idx]
    λi = λ[i,idx]

    scatter!(sp1, Vi, imag.(λi),
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

for i = 1:size(λ, 1)

    idx = findall(x -> abs(x) < 500, λ[i,:])
    Vi = Vinf[idx]
    λi = λ[i,idx]

    scatter!(sp2, Vi,
        # real.(λi),
        real.(λi)./abs.(λi)*100,
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800), show=true)
