# functions to plan in given environment 
using Graphs
mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # path length in euclidean space
end
function plan(world::PlanningProblem)::Path
    # empty undirected graph
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # Stores 3D state space points corresponding to graph vertices
    controls = Vector{Vector{Float64}}()  # Stores control inputs corresponding to edges
    # known initial point
    init_point = world.env.init
    add_vertex!(g)
    push!(vertices, init_point)
    # create implicit reward function around ROIs
        # function reward =  
    # Set up environment
    goal = world.env.goal
    lim = world.env.lim  # Bounds
    obstacles = world.env.obstacles
    max_iter = 10000    # max iter
    dt = 0.1            
    # Runtime loop
    for i in 1:max_iter
        # Sample a random point within bounds
        rand_point = if rand() > 0.5
            goal
        else
            samplePoint(lim)            
        end
        # Find the nearest vertex in the graph
        nearest_idx = argmin([norm(rand_point - v[1:3]) for v in vertices])  # Use only position
        nearest_x = vertices[nearest_idx]
        #nearest_point = nearest_x[1:3]

        # Find best control towards the random point
        ustar = getBestControl(nearest_x, rand_point, dt, world.agents[1], world.env.windField)

        # Propagate dynamics to generate a new point
        new_x = propagate(nearest_x, ustar, dt, world.agents[1].EOM, world.env.windField)
        new_point = new_x[1:3]
        # Check if the new point is valid
        if isPointValid(new_point, lim)
            add_vertex!(g)
            push!(vertices, new_x)
            push!(controls, ustar)  # Save control corresponding to this edge
            add_edge!(g, nearest_idx, length(vertices))

            # Check if the goal is reached
            if norm(new_point - goal) <= 1.0
                return extract_path(g, vertices, controls, init_point, goal)
            end
        end
    end

    error("RRT failed to find a path within the maximum number of iterations")
end
function extract_path(
    g::SimpleGraph,
    vertices::Vector{Vector{Float64}},
    controls::Vector{Vector{Float64}},
    start::Vector{Float64},
    goal::Vector{Float64}
)::Path
    # Extract path from graph using backtracking
    wp = [goal]
    control_seq = Vector{Vector{Float64}}()
    current_idx = length(vertices)

    while vertices[current_idx] != start
        prev_idx = outneighbors(g, current_idx)[1]
        push!(wp, vertices[prev_idx])
        push!(control_seq, controls[current_idx-1])  # Controls are indexed by edge
        current_idx = prev_idx
    end
    # Reverse the waypoints and controls to maintain correct order
    wp = reverse(wp)
    control_seq = reverse(control_seq)

    # Compute path length
    path_length = sum(norm(wp[i+1] - wp[i]) for i in 1:(length(wp)-1))

    return Path(wp, control_seq, path_length)
end
# get a random point within bounds
function samplePoint(l::Bounds)::Vector{Float64}
    x = rfinbounds(l.xmin, l.xmax)
    y = rfinbounds(l.ymin, l.ymax)
    z = rfinbounds(l.zmin, l.zmax)
    return [x, y, z]
end
# get a random control within bounds
function randControl(ubounds::Vector{Vector{Float64}})::Vector{Float64}
    rf = zeros(Float64, length(ubounds))
    for i in 1:length(ubounds)
        rf[i] = rfinbounds(ubounds[i][1], ubounds[i][2])
    end
    return rf
end
# try a bunch of controls to see what works
function getBestControl(xk::Vector{Float64}, xkp1::Vector{Float64}, dt::Float64, agent::Agent, wind::Function)::Vector{Float64}
    best_control = randControl(agent.controlBounds)
    best_score = -Inf
    for i in 1:10
        u = randControl(agent.controlBounds)
        new_point = propagate(xk,u,dt,agent.EOM,wind)
        score = scoreControl(new_point, xkp1)
        if score > best_score
            best_score = score
            best_control = u
        end
    end
    return best_control
end
# score a control 
function scoreControl(new_point::Vector{Float64}, goal_point::Vector{Float64})::Float64
    return -norm(new_point - goal_point)  # closer is better
end
# check if point is valid (airspeed threshold)
function isPointValid(point::Vector{Float64}, _)::Bool
    return true
end
# propagate system using RK4 dynamics
function propagate(x::Vector{Float64},u::Vector{Float64},dt::Float64,eom::Function,wind::Function)::Vector{Float64}
    xkp1 = RK4_int_wctrl(dt,x,u,eom,wind)
    return xkp1
end
# Runge-Kutta 4th-order integrator
function RK4_int_wctrl(dt, state, u, fdot::Function, wind::Function)::Vector{Float64}
    # assume ZOH
    k1 = dt * fdot(state,u,wind)
    k2 = dt * fdot(state + 0.5 * k1, u,wind)
    k3 = dt * fdot(state + 0.5 * k2, u,wind)
    k4 = dt * fdot(state + k3,u,wind)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end