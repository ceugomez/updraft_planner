include("./src/generic.jl");  # environment and problem definitions from updraft_planner
include("./src/planner.jl");  # planner types from updraft_planner
include("./src/winds/generic.jl") # generic methods (wind field) from 12dofAC
include("./src/optim/optimizer.jl") # optimizer and info field from InfoOptimizer
include("./src/dynamics.jl")    # dynamics functions
include("./src/visualizer.jl") #visualization functions
function windField_planner(pos::Vector{Float64})
    x,y,z=pos[1:3]
    return [0.0,0.0,0.0]
end
function simulate_dynamics(state::Vector{Float64}, dt::Float64, interval::Float64, eom::Function, wind::Function)::Path
    # Initialize variables
    t = 0.0  # Start time
    trajectory = [state]  # Store the trajectory states

    # Define zero control input
    zero_control = [-0.01, 0.00]  # Assuming the control vector has (ψ_dot, γ_dot)

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
        x0_AC = [100.0, 0.0, 4500.0, 25.0, -0.01, -0.056];                                  # initial point [x y z Vi ψ γ ]
        xg_AC = [500.0, 0.0, 125.0, 20.0, 0.0, 0.0]                             # goal point [x y z Vi _ _]
        controlLim_AC = [[-0.01, 0.01],[-0.05, -0.6]]   # [ψdot, γdot]
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
        plot_dubins_glider_path(path_AC, 1.0)

    # simulate simply to check dynamics
        # hist = simulate_dynamics(x0_AC, 1.0, 1000.0, dubins_glider_dynamics, windField_planner)    
        # plot_dubins_glider_path(hist, 1.0)        

end
