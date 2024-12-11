
function simpleIntegrator(x,u,wind)
    # simple integrator, placeholder
    return x+0.1*u
end
function dubins_glider_dynamics(state::Vector{Float64}, control::Vector{Float64}, dt::Float64)::Vector{Float64}
    # Unpack the state variables
    x, y, z, θ, γ, v = state

    # Unpack the control inputs
    ω, ϕ = control

    # Physical constants
    g = 9.81  # Gravity acceleration (m/s^2)
    glide_ratio = 15.0  # Glide ratio (lift-to-drag)

    # Compute lift and drag forces
    drag = g / sqrt(1 + glide_ratio^2)  # Drag component
    lift = drag * glide_ratio           # Lift component

    # Update velocity considering drag
    v_new = v - drag * dt

    # Ensure airspeed is positive
    v_new = max(v_new, 0.1)

    # Compute the next state using glider kinematics
    dx = v_new * cos(θ) * cos(γ) * dt  # Horizontal x component
    dy = v_new * sin(θ) * cos(γ) * dt  # Horizontal y component
    dz = v_new * sin(γ) * dt           # Vertical z component
    dθ = ω * dt                        # Change in heading angle
    dγ = (lift / v_new - g / v_new) * dt + ϕ * dt  # Change in flight path angle

    # Update the state
    next_state = [
        x + dx,
        y + dy,
        z + dz,
        θ + dθ,
        γ + dγ,
        v_new
    ]

    # Wrap heading angle θ to [-π, π]
    next_state[4] = mod(next_state[4] + π, 2π) - π

    # Clamp γ (flight path angle) to a physical range (e.g., -π/6 to π/6 for this model)
    γ_min, γ_max = -π/6, π/6
    next_state[5] = clamp(next_state[5], γ_min, γ_max)

    return next_state
end