
function simpleIntegrator(x,u,wind)
    # simple integrator, placeholder
    return u
end
function dubins_glider_dynamics(state::Vector{Float64}, u::Vector{Float64}, wind::Function)::Vector{Float64}
    # Unpack the state variables
    x, y, z, v, ψ, γ,  = state  # x, y, z: positions; v: velocity; ψ: heading angle; γ: glide slope angle; 

    # Unpack the control inputs
    ω, ϕ = u  # ω: change in heading angle (rate), ϕ: change in glide slope angle (rate)

    Vw = wind([x,y,z])[3];
    # Physical constants
    g = 9.81        # grav. accel [m/s^2]
    S = 0.628       # wing  area [m^2]
    m = 5.74        #  mass [kg]
    Cd0 = 0.0140    # parasite drag coefficient
    Clα = 6.196683  # lift-curve slope [1/rad]
    ρ = 1.025       # air density (assume const) [kg/m^3]
    
    # aero coeffs
    α = ϕ       # Assume glide slope angle as proxy for angle of attack
    Cl = Clα * α  # Lift coefficient (linear approximation)
    Cd = Cd0 + (Cl^2 / (π * 4.3 * 0.8))  # Drag coefficient (polar model; e ~ 0.8, AR ~ 4.3)
    D = 0.5 * ρ * v^2 * S * Cd  # Drag force
    #L = 0.5 * ρ * v^2 * S * Cl  # Lift force

    # get state derivatives
    dx = v * cos(γ) * cos(ψ)       # xdot
    dy = v * cos(γ) * sin(ψ)       # ydot
    dz = v * sin(γ) + Vw          # zdot
    dv = -g * sin(γ) - D / m       # vdot
    dψ = ω                         # Change in heading angle
    dγ = ϕ                         # Change in glide slope angle
    # return state derivative
    return [dx, dy, dz, dv, dψ, dγ]
end
