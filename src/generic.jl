using Distributions, Parameters, Atmosphere 
# problem boundaries in 3-space
struct Bounds
    xmin::Float16       # distance, m
    xmax::Float16       # hopefully these are self-explanatory
    ymin::Float16
    ymax::Float16
    zmin::Float16
    zmax::Float64
end

#environment definition, with obstacles, goals, and  
mutable struct Environment
    lim::Bounds
    obstacles::Vector{Vector{Float64}}  # obstacles, vec[ vec[vertices] ]
    init::Vector{Float64}               # Initial Point,     [x,y,z]
    ROI::Vector{Vector{Float64}}        # Regions to sample, [ [x,y,z] ] 
    goal::Vector{Float64}               # Finish point       [x,y,z]
    windField::Function
end
# agent definition, with equations of motion and physical dimensions
struct Agent
    EOM::Function
    controlBounds::Vector{Vector{Float64}}
    wsBounds::Bounds
end
# structure to store problem information
mutable struct PlanningProblem
    env::Environment             # environent to plan in
    agents::Vector{Agent}         # agents to use 
end

function windField_planner(pos::Vector{Float64})::Vector{Float64}
    x, y, z = pos[1:3]
    
    # Define parameters for the wind field regions
    region1_center = [975.0, 440.0, 50.0]
    region1_radius = 750.0
    region1_strength = -1.0  # m/s

    region2_center = [-335, -75.5, 50.0]
    region2_radius = 450.0
    region2_strength = 2.0  

    region3_center = [1500.0, 1500.0, 50.0]
    region3_radius = 150.0
    region3_strength = 3.0  

    # Compute distances to the centers of the regions
    dist1 = norm([x, y, z] - region1_center)
    dist2 = norm([x, y, z] - region2_center)
    dist3 = norm([x, y, z] - region3_center)

    # Compute contributions from each region (use Gaussian-like decay)
    region1_flow = region1_strength * exp(-dist1^2 / (2 * region1_radius^2))
    region2_flow = region2_strength * exp(-dist2^2 / (2 * region2_radius^2))
    region3_flow = region3_strength * exp(-dist3^2 / (2 * region3_radius^2))

    # Superposition of contributions
    wz_raw = region1_flow + region2_flow + region3_flow

    # Scale updrafts logarithmically with altitude, capped at maximum value
    z_min = 1.0  # Prevent log(0)
    z_max = 500.0  # Adjust based on the maximum altitude in your domain
    log_scaling = log(max(z, z_min)) / log(z_max)  # Scales between 0 and 1
    wz = wz_raw * clamp(log_scaling, 0.0, 1.0)  # Clamp to avoid overshooting

    # Optional: Add small horizontal wind flow
    wx = 0.0  # Modify as needed
    wy = 0.0  # Modify as needed

    return [wx, wy, wz]
end

function simulate_dynamics(state::Vector{Float64}, dt::Float64, interval::Float64, eom::Function, wind::Function)::Path
    # Initialize variables
    t = 0.0  # Start time
    trajectory = [state]  # Store the trajectory states

    # Define zero control input
    zero_control = [-0.001, 0.00]  # Assuming the control vector has (ψ_dot, γ_dot)

    # Simulate until the interval limit
    while t < interval
        # Get the current state
        current_state = trajectory[end]

        # Propagate dynamics with zero control
        next_state = propagate(current_state, zero_control, dt, eom, wind)

        # Update trajectory
        push!(trajectory, next_state)

        # Increment time
        t += dt
    end
    temp = Path(trajectory,[[0.0]],0.0)
    return temp
end