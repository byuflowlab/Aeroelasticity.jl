"""
    Section

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct Section end

"""
    Section()

Initialize a model of type [`Section`](@ref)
"""
Section()

(::Section)(resid, dx, x, y, p, t) = section_residual!(resid, dx, x, y, p, t)

number_of_states(::Section) = 4

number_of_parameters(::Section) = 5

# --- Internal Methods --- #

# residual function
function section_residual!(resid, dx, x, y, p, t)

    dh, dθ, dhdot, dθdot = dx
    h, θ, hdot, θdot = x
    L, M = y
    kh, kθ, m, Sθ, Iθ = p   

    resid[1] = dh - hdot
    resid[2] = dθ - θdot
    resid[3] = m*dhdot + Sθ*dθdot + kh*h + L
    resid[4] = Sθ*dhdot + Iθ*dθdot + kθ*θ - M

    return resid
end

# airfoil local linear/angular velocities
function section_velocities(U, θ, hdot, θdot)
    u = U
    v = U*θ + hdot
    ω = θdot
    return SVector(u, v, ω)
end

# airfoil local linear/angular accelerations
function section_accelerations(dhdot, dθdot)
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    return SVector(udot, vdot, ωdot)
end

# airfoil coordinates (for plotting)
function section_coordinates(h, θ; 
    a=0, 
    b=0.5,
    xcoord = [1.0, 0.993844, 0.975528, 0.945503, 0.904508, 0.853553, 0.793893, 0.726995, 
        0.654508, 0.578217, 0.5, 0.421783, 0.345492, 0.273005, 0.206107, 0.146447, 
        0.095492, 0.054497, 0.024472, 0.006156, 0.0, 0.006156, 0.024472, 0.054497, 
        0.095492, 0.146447, 0.206107, 0.273005, 0.345492, 0.421783, 0.5, 0.578217, 
        0.654508, 0.726995, 0.793893, 0.853553, 0.904508, 0.945503, 0.975528, 0.993844, 1.0],
    ycoord = [0.00126, 0.00212, 0.004642, 0.008658, 0.013914, 0.020107, 0.026905, 0.033962, 
        0.040917, 0.047383, 0.05294, 0.057148, 0.059575, 0.059848, 0.057714, 0.053083, 
        0.046049, 0.036867, 0.025893, 0.013503, 0.0, -0.013503, -0.025893, -0.036867, 
        -0.046049, -0.053083, -0.057714, -0.059848, -0.059575, -0.057148, -0.05294, 
        -0.047383, -0.040917, -0.033962, -0.026905, -0.020107, -0.013914, -0.008658, 
        -0.004642, -0.00212, -0.00126],
    )

    xplot = similar(xcoord)
    yplot = similar(ycoord)
    for i = 1:length(xcoord)
        xplot[i] = (xcoord[i] - 0.5 - a)*2*b*cos(θ) - ycoord[i]*2*b*sin(θ)
        yplot[i] = (xcoord[i] - 0.5 - a)*2*b*sin(θ) + ycoord[i]*2*b*cos(θ) + h
    end

    return xplot, yplot
end