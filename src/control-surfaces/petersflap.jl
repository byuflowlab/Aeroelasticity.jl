"""
    PetersFlap <: NoStateModel

Linear, steady-state, two-dimensional control surface model with parameters
``c_{n,\\delta}``, ``c_{a,\\delta}``, and ``c_{m,\\delta}``.
"""
struct PetersFlap <: NoStateModel end

"""
    PetersFlap()

Initialize an object of type [`PetersFlap`](@ref)
"""
PetersFlap()

# --- Traits --- #

number_of_parameters(::Type{PetersFlap}) = 3

# --- Convenience Methods --- #

function set_parameters!(p, model::PetersFlap; cnd, cad, cmd)

    d
    ε
    δ

    p[1] = cnd
    p[2] = cad
    p[3] = cmd

    return p
end

function separate_parameters(model::PetersFlap, p)

    return (cnd=p[1], cad=p[2], cmd=p[3])
end

# --- Internal Methods --- #

function petersflap_loads(b, U, ρ, cnδ, caδ, cmδ, δ)

    T = petersflap_angular_deflection_coefficients(θ, n)
    h = T*δ
    Σ = (1:N)'*h
    L = 2*pi*ρ*U^2*Σ
    D =
    M = 2*pi*ρ*U^2*h[1]

    n = 2
    hn = petersflap_displacement_field
    2*sin(n*θ)/(n*pi) * ε + b/pi*((sin(n-1)*θ)/(n-1) - (1+d)*sin(n*θ)/n + sin((n+1)*θ)/(n+1))*δ

    L = ρ*U^2*b*cnδ*δ
    M = ρ*U^2*b^2*cmδ*δ
    A/(2*pi*ρ*b) = -(v + d*ω - λ0 - u*α0)^2 - (v + d*ω - λ0 - u*α0)*u/b*Σ +
        ω1*U2*h[1] -1/2*vdot*h[1] + b/4*ωdot*h[2]

    A/(2*pi*ρ*b) = -(v + d*ω - λ0 - u*α0)^2 + (Σ/b*u*v + (Σ/b*d + h[1])*u*ω - Σ/b*u*λ0 - u^2*α0*Σ/b - 1/2*vdot*h[1]) +
         + b/4*ωdot*h[2]


    return L, D, M
end

function petersflap_displacement_field(b, θ, ε, δ, n)
    if n == 0
        return θ/pi*ε + b/pi*(sin(θ) - (1+d)/2*θ)*δ
    elseif n == 1
        return 2*sin(θ)/π*ε + b/pi*(θ - (1+d)*sin(θ) + sin(2*θ)/2)*δ
    else
        return 2*sin(n*θ)/(n*pi) * ε + b/pi*((sin(n-1)*θ)/(n-1) -
            (1+d)*sin(n*θ)/n + sin((n+1)*θ)/(n+1))*δ
    end
end



function petersflap_linear_deflection_coefficients(θ, n)
    return SVector(ntuple(n -> petersflap_linear_deflection_coefficient(θ, n), n)...)
end

function petersflap_linear_deflection_coefficient(θ, n)
    if iszero(n)
        return θ/pi
    else
        return 2*sin(θ)/(n*pi)
    end
end

function petersflap_angular_deflection_coefficients(θ, n)
    return SVector(ntuple(n -> petersflap_angular_deflection_coefficient(θ, n), n)...)
end

function petersflap_angular_deflection_coefficient(θ, n)
    if iszero(n)
        return sin(θ) - θ*cos(θ)
    else
        return b/pi*sin((n+1)*θ)/(n+1) - 2/n*cos(θ)*sin(n*θ) + sin((n-1)*θ)/(n-1)
    end
end
