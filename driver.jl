include("./src/generic.jl");  # environment and problem definitions from updraft_planner
include("./src/planner.jl");  # planner types from updraft_planner
include("./src/winds/generic.jl") # generic methods (wind field) from 12dofAC
include("./src/optim/optimizer.jl") # optimizer and info field from InfoOptimizer
include("./src/dynamics.jl")    # dynamics functions
include("./src/visualizer.jl") #visualization functions
function windField_planner(pos::Vector{Float64})::Vector{Float64}
    x, y, z = pos[1:3]
    
    # Define parameters for the wind field regions
    region1_center = [50.0, 50.0, 50.0]
    region1_radius = 40.0
    region1_strength = 2.0  # Maximum upward flow

    region2_center = [-50.0, -50.0, 50.0]
    region2_radius = 40.0
    region2_strength = -2.0  # Maximum downward flow

    # Compute distances to the centers of the regions
    dist1 = norm([x, y, z] - region1_center)
    dist2 = norm([x, y, z] - region2_center)

    # Compute contributions from each region (use Gaussian-like decay)
    region1_flow = region1_strength * exp(-dist1^2 / (2 * region1_radius^2))
    region2_flow = region2_strength * exp(-dist2^2 / (2 * region2_radius^2))

    # Combine contributions
    wz = region1_flow + region2_flow

    # Add a small horizontal wind flow (optional)
    wx = 0.5 * sin(x / 20.0)  # Horizontal flow component along x
    wy = 0.5 * cos(y / 20.0)  # Horizontal flow component along y

    return [wx, wy, wz]
end
function simulate_dynamics(state::Vector{Float64}, dt::Float64, interval::Float64, eom::Function, wind::Function)::Path
    # Initialize variables
    t = 0.0  # Start time
    trajectory = [state]  # Store the trajectory states

    # Define zero control input
    zero_control = [0.0, 0.00]  # Assuming the control vector has (ψ_dot, γ_dot)

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
begin
    # define environment - 3DOF AC model 
        wslim_AC = Bounds(-2000.0,2000.0,-2000.0,2000.0,0.0,5000.0);
        x0_AC = [-1000.0, 0.0, 4500.0, 25.0, 0.0, -0.056];                                  # initial point [x y z Vi ψ γ ]
        xg_AC = [-550, 450.0, 1500.0, 0.0, 0.0, 0.0]                             # goal point [x y z Vi _ _]
        controlLim_AC = [[-0.3, 0.3],[-0.06, 0.06]]   # [ψdot, γdot]
        env_AC = Environment(wslim_AC, [[]], x0_AC, [[]],xg_AC , windField_planner)    # ws dimensions, obstacles, init, roi (added iteratively), goal
        agents_AC = [Agent(dubins_glider_dynamics, controlLim_AC, wslim_AC)];                               # agent definition 
        planprob_AC = PlanningProblem(env_AC,agents_AC)

    # define environment - Simple Integrator
        # wslim_SI = Bounds(-30.0,30.0,-30.0,30.0,0.0,2.0);
        # x0_SI = [28.0, 0.0, 1.9]
        # xg_SI = [0.0, 0.0, 0.0]
        # controlLim_SI = [[-1.0,1.0],[-1.0, 1.0], [-1.0, 1.0]]
        # env_SI = Environment(wslim, [[]], x0_SI, [[]],xg_SI , WindField) # ws dimensions, obstacles, init, roi (added iteratively), goal
        # agents_SI = [Agent(simpleIntegrator, controlLim_SI, wslim)];                               # agent definition 
        # planprob_SI = PlanningProblem(env_SI, agents_SI)                                                     # 

    # define and solve ROI optimization
        # OptProb = OptimProblemStatement(length(planprob.agents),         # number of agents to optimize for 
        #                         [wslim.xmin wslim.xmax; wslim.ymin wslim.ymax; wslim.zmin wslim.zmax], 
        #                         "Cross-Entropy",           # optimizer method
        #                         Vector{Vector{Float64}}(), # vector of positions, empty
        #                         5.0                        # exclusion radius for each ROI 
        #                         ); 
    # define and solve planning problem
        #planprob.env.ROI = OptimizeEnvironment(OptProb, 1e-2);
        path_AC = plan(planprob_AC);   
        # plotPath3d(path_AC, planprob_AC)
        # plot_dubins_glider_path(path_AC, 1.0)

    # simulate simply to check dynamics
        #  hist = simulate_dynamics(x0_AC, 1.0, 1000.0, dubins_glider_dynamics, windField_planner)    
        #  plot_dubins_glider_path(hist, 1.0)        

end
