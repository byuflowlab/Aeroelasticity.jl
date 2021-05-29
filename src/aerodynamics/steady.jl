struct Steady <: NoStateModel end

# --- Model Traits --- #
number_of_parameters(::Type{Steady}) = 3

# --- Coupled Model Traits --- #
inplaceness(::Type{Steady}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{Steady}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{Steady}, ::Type{TypicalSection}) = Varying()

# --- Coupled Model Methods --- #
function get_inputs(aero::Steady, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract parameters
    b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # calculate aerodynamic loads
    L, M = steady_loads(b, U, ρ, θ)
    # return inputs
    return SVector(L, M)
end

function get_input_state_jacobian(aero::Steady, stru::TypicalSection, u, p, t)
    # extract parameters
    b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0; 2*pi*ρ*b*U^2 0; 0 0; 0 0]
end

# TODO: Parameter jacobian

# --- Internal --- #
function steady_loads(b, U, ρ, θ)
    L = 2*pi*ρ*b*U^2*θ
    M = 0 # quarter chord moment
    return SVector(L, M)
end
