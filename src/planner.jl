# functions to plan in given environment 
using Graphs
mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # path length in euclidean space
end
function plan(world::PlanningProblem)::Vector{Vector{Float64}}
    # Initialize an empty undirected graph
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # Stores 3D points corresponding to graph vertices
    
    # Add initial point
    init_point = world.env.init
    add_vertex!(g)
    push!(vertices, init_point)
    
    # Set up environment
    goal = world.env.goal
    lim = world.env.lim  # Bounds
    obstacles = world.env.obstacles
    max_iter = 10000  # Maximum number of iterations
    step_size = 1.0  # Step size
    dt = 0.1         # Time resolution for propagation
    
    # Runtime loop
    for i in 1:max_iter
        # Sample a random point within bounds
        if (rand()>0.5)
            rand_point = goal
        else
            rand_point = samplePoint(lim)
        end
        
        # Find the nearest vertex in the graph
        nearest_idx = argmin([norm(rand_point - v) for v in vertices])
        nearest_point = vertices[nearest_idx]
        
        # Find best control towards the random point
        ustar = getBestControl(nearest_point, rand_point, dt, world.agents[1], world.env.windField)  # Assuming no wind for now
        
        # Propagate dynamics to generate a new point
        new_point = propagate(nearest_point, ustar, dt, world.agents[1].EOM, world.env.windField)
        
        # Check if the new point is valid
        if isPointValid(new_point, lim, obstacles)
            add_vertex!(g)
            push!(vertices, new_point)
            add_edge!(g, nearest_idx, length(vertices))
            
            # Check if the goal is reached
            if norm(new_point - goal) <= step_size
                return extract_path(g, vertices, init_point, goal)
            end
        end
    end
    
    error("RRT failed to find a path within the maximum number of iterations")
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

function samplePoint(l::Bounds)::Vector{Float64}
    x = rfinbounds(l.xmin, l.xmax)
    y = rfinbounds(l.ymin, l.ymax)
    z = rfinbounds(l.zmin, l.zmax)
    return [x, y, z]
end

function randControl(ubounds::Vector{Vector{Float64}})::Vector{Float64}
    rf = zeros(Float64, length(ubounds))
    for i in 1:length(ubounds)
        rf[i] = rfinbounds(ubounds[i][1], ubounds[i][2])
    end
    return rf
end

function getBestControl(xk::Vector{Float64}, xkp1::Vector{Float64}, dt::Float64, agent::Agent, wind::Function)::Vector{Float64}
    best_control = randControl([[0.0, 1.0], [0.0, 1.0]]) 
    best_score = -Inf
    for i in 1:100
        u = randControl(agent.controlBounds)
        new_point = propagate(xk,u,dt,agents[1].EOM, wind)
        score = scoreControl(new_point, xkp1)
        if score > best_score
            best_score = score
            best_control = u
        end
    end
    return best_control
end

function scoreControl(new_point::Vector{Float64}, goal_point::Vector{Float64})::Float64
    # Score based on proximity to the goal and other objective
    return -norm(new_point - goal_point)  # closer is better
end

function isPointValid(point::Vector{Float64}, lim::Bounds, obstacles::Vector{Vector{Float64}})::Bool
    # Check if the point is within bounds
    if !(lim.xmin <= point[1] <= lim.xmax &&
         lim.ymin <= point[2] <= lim.ymax &&
         lim.zmin <= point[3] <= lim.zmax)
        return false
    end
    
    # Check for collisions with obstacles
    for obstacle in obstacles
        if in_obstacle(point, obstacle)
            return false
        end
    end
    
    return true
end

function in_obstacle(point::Vector{Float64},_)::Bool
    # Example assumes obstacles are axis-aligned bounding boxes
    return false
end

function propagate(start::Vector{Float64}, u::Vector{Float64}, dt::Float64, eom::Function, wind::Function)::Vector{Float64}
    # Propagate dynamics from the start point with control `u` over time `dt`
    return eom(start, u,wind)  # Assuming EOM function integrates for one timestep
end
