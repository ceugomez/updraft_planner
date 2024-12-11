include("./src/generic.jl");  # environment and problem definitions from updraft_planner
include("./src/planner.jl");  # planner types from updraft_planner
include("./src/winds/generic.jl") # generic methods (wind field) from 12dofAC
include("./src/optim/optimizer.jl") # optimizer and info field from InfoOptimizer
include("./src/dynamics.jl")    # dynamics functions
include("./src/aircraft_sim/aircraft_eom.jl")
include("./src/visualizer.jl") #visualization functions
begin
    # define environment 
    wslim = Bounds(-30.0,30.0,-30.0,30.0,0.0,2.0);
    env = Environment(wslim, [[]], [28.0, 0.0, 1.9], [[]], [0.0, 0.0, 0.0], WindField) # ws dimensions, obstacles, init, roi (added iteratively), goal
    controlLim = [[-0.05, 0.05],[-0.05, 0.05],[-0.05, 0.05],[0,1]] # [δe, δa, δr, δt]
    
    agents = [Agent(acEOM, controlLim, wslim)];                               # agent definition (point for now)
    planprob = PlanningProblem(env, agents)    # 
    # define and solve ROI optimization
        OptProb = OptimProblemStatement(length(planprob.agents),         # number of agents to optimize for 
                                [wslim.xmin wslim.xmax; wslim.ymin wslim.ymax; wslim.zmin wslim.zmax], 
                                "Cross-Entropy",           # optimizer method
                                Vector{Vector{Float64}}(), # vector of positions, empty
                                5.0                        # exclusion radius for each ROI 
                                ); 
    # define and solve planning problem
        #planprob.env.ROI = OptimizeEnvironment(OptProb, 1e-2);
        path = plan(planprob);   
        #print(typeof(path))       
        plotPath3d(path, planprob)

end