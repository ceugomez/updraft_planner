include("generic.jl");  # environment and problem definitions from updraft_planner
include("planner.jl");  # planner types from updraft_planner
include("12dofAC/generic.jl") # generic methods (wind field) from 12dofAC
include("esa_InfoOptimizer/optimizer.jl") # optimizer and info field from InfoOptimizer

begin
    # define environment 
    wslim = Bounds(-30.0,30.0,-30.0,30.0,0.0,2.0);
    env = Environment(wslim, [[]], [0.0, 0.0, 1.9], [[]], [0.0, 0.0, 0.0] ) # ws dimensions, obstacles, init, roi (added iteratively), goal
    agents = [Agent(acEOM)];                               # agent definition (point for now)
    planprob = PlanningProblemDefn(env, agents)    # 
    # define and solve ROI optimization
        OptProb = OptimProblemStatement(length(planprob.agents),         # number of agents to optimize for 
                                [wslim.xmin wslim.xmax; wslim.ymin wslim.ymax; wslim.zmin wslim.zmax], 
                                "Cross-Entropy",           # optimizer method
                                Vector{Vector{Float64}}(), # vector of positions, empty
                                5.0                        # exclusion radius for each ROI 
                                ); 
    planprob.env.ROI = OptimizeEnvironment(OptProb, 1e-2)
    println(planprob.env.ROI)
    plotPath2d(planprob.env.ROI)
          

end