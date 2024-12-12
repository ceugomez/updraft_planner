include("./src/generic.jl");  # environment and problem definitions from updraft_planner
include("./src/planner.jl");  # planner types from updraft_planner
include("./src/winds/generic.jl") # generic methods (wind field) from 12dofAC
include("./src/optim/optimizer.jl") # optimizer and info field from InfoOptimizer
include("./src/dynamics.jl")    # dynamics functions
include("./src/visualizer.jl") #visualization functions
begin
    # define environment - 3DOF AC model 
        wslim_AC = Bounds(-300.0,300.0,-300.0,300.0,0.0,200.0);
        x0_AC = [28.0, 0.0, 1.9, 18.0, 0.0, 0.0];                                  # initial point [x y z Vi ψ γ ]
        xg_AC = [0.0, 0.0, 0.0, 15.0, 0.0, 0.0]                             # goal point [x y z Vi _ _]
        controlLim_AC = [[-0.05, 0.05],[-0.05, 0.05]]   # [ψdot, γdot]
        env_AC = Environment(wslim_AC, [[]], x0_AC, [[]],xg_AC , WindField)    # ws dimensions, obstacles, init, roi (added iteratively), goal
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
        plotPath3d(path_AC, planprob_AC)

end