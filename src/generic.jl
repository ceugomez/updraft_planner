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

