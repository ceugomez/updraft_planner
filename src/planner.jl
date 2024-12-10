# functions to plan in given environment 
using Graphs
mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # path length in euclidean space
end

function plan(world::PlanningProblem)::Path
    # initialize empty undirected graph
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # Stores 3D points corresponding to graph vertices
    # add initial point
    init_point = world.env.init
    print(typeof(init_point))
    add_vertex!(g)
    push!(vertices, init_point)
    # set up environment     
    goal = world.env.goal
    lim = world.env.lim  # Bounds 
    obstacles = world.env.obstacles
    max_iter = 1000  # max no. of iterations
    dt = 1.0  # timestep 
    # runtime loop   
    for i in 1:max_iter
        # sample a random point within bounds
        rand_point = samplePoint(world.env.lim)
        # find the nearest vertex in the graph
        nearest_idx = argmin([norm(rand_point - v) for v in vertices])
        nearest_point = vertices[nearest_idx]
        # propagate towards the random point using known dynamics
            ustar = findBestControl(nearest_point, sample, )
            new_point = world.agents[1].EOM(nearest_point, rand_point, step_size, world.agents[1].EOM)
        # find best control
        add_vertex!(g)
        push!(vertices, new_point)
        add_edge!(g, nearest_idx, length(vertices))
        # Check if goal is reached
        if norm(new_point - goal) <= step_size
            return extract_path(g, vertices, init_point, new_point)
        end

    end

    error("RRT failed")
end

function extract_path(g::SimpleGraph, vertices::Vector{Vector{Float64}}, start::Vector{Float64}, goal::Vector{Float64})::Vector{Vector{Float64}}
    # Extract path from graph using backtracking
    path = [goal]
    current_idx = length(vertices)
    while vertices[current_idx] != start
        current_idx = outneighbors(g, current_idx)[1]
        push!(path, vertices[current_idx])
    end
    return reverse(path)
end
# sample random point within domain
function samplePoint(l::Bounds)::Vector{Float64}
    x = rfinbounds(l.xmin, l.xmax)
    y = rfinbounds(l.ymin, l.ymax)
    z = rfinbounds(l.zmin, l.zmax)
    return [x, y, z]
end
# get random control within control bounds
function randControl(ubounds::Vector{Vector{Float64}})
    rf = zeros(1, length(ubounds))
    for i in i:length(ubounds)
        rf[i] = rfinbounds(ubounds[i][1], ubounds[i][2])
    end
    return rf
end
# get the best control in the direction of xkp1 from x
function getBestControl(xk::Vector{Float64}, xkp1::Vector{Float64at}, agent::Agent, wind::Function)
    u = Vector{Float64}(0.0,0.0)
    for i in 1:100
        u = randControl
        x = agent.EOM(xk )
    end
end
function scoreControl()
    # score based on potential energy gain + progress towards objectives 
    # Bellman–Ford algorithm
    
end
function isPointValid()
    # check if point is in bounds
end