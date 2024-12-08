using Distributions, Parameters, Atmosphere 

# structure to store problem information
@with_kw mutable struct ProblemDfn
    lim::Bounds
    env::Environment
    agent::Agent
end
# problem boundaries in 3-space
struct Bounds
    xmin::Float16       # distance, m
    xmax::Float16       # hopefully these are self-explanatory
    ymin::Float16
    ymax::Float16
    zmin::Float16
    zmax::Float64
end
# environment definition, with obstacles, goals, and  
struct Environment
    obstacles::Vector{Vector{Float64}}  # obstacles, vec[ vec[vertices] ]
    init::Vector{Float64}               # Initial Point,     [x,y,z]
    ROI::Vector{Vector{Float64}}        # Regions to sample, [ [x,y,z] ] 
    goal::Vector{Float64}               # Finish point       [x,y,z]
end
# agent definition, with equations of motion and physical dimensions
struct Agent
    EOM::Function


end

