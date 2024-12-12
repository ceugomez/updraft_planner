using Graphs
using NearestNeighbors
using Distributed

mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # Vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # Path length in Euclidean space
end
function plan(world::PlanningProblem)::Path
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # Stores state space points corresponding to graph vertices
    controls = Vector{Vector{Float64}}()  # Stores control inputs corresponding to edges
    parents = Dict{Int, Int}()            # Store parent indices for path extraction

    init_point = world.env.init
    add_vertex!(g)
    push!(vertices, init_point)

    goal = world.env.goal
    lim = world.env.lim
    max_iter = 1000
    dt = 1.0
    kdtree_update_frequency = 100  # Update KDTree every 100 iterations
    kdtree = KDTree(hcat([init_point[1:3]]...))  # Initial KDTree

    for i in 1:max_iter
        # Update KDTree periodically
        if i % kdtree_update_frequency == 0
            kdtree = KDTree(hcat([v[1:3] for v in vertices]...))
        end

        # Sample a random point
        goal_sampling_prob = min(0.5, length(vertices) / max_iter)  # close towards goal sampling over time
        rand_point = if rand() < goal_sampling_prob
            goal[1:3]
        else
            samplePoint(lim)
        end

        # Find the nearest vertex using KDTree
        nearest_idx = knn(kdtree, rand_point, 1)[1][1]
        nearest_x = vertices[nearest_idx]

        # Find best control towards the random point
        ustar = getBestControl(nearest_x, rand_point, dt, world.agents[1], world.env.windField)

        # Propagate dynamics to generate a new point
        new_x = propagate(nearest_x, ustar, dt, world.agents[1].EOM, world.env.windField)

        # Check if the new point is valid
        if isPointValid(new_x, lim)
            add_vertex!(g)
            push!(vertices, new_x)
            push!(controls, ustar)
            add_edge!(g, nearest_idx, length(vertices))
            parents[length(vertices)] = nearest_idx  # Update parent index

            # Check if the goal is reached
            if norm(new_x[1:3] - goal[1:3]) <= 10.0
                println("Reached goal!")
                return extract_path(g, vertices, controls, parents, init_point, goal, true)
            end
        end
    end

    println("RRT failed to find a path within the maximum number of iterations")
    return extract_path(g, vertices, controls, parents, init_point, goal, false)
end

function extract_path(
    g::SimpleGraph,
    vertices::Vector{Vector{Float64}},
    controls::Vector{Vector{Float64}},
    parents::Dict{Int, Int},
    start::Vector{Float64},
    goal::Vector{Float64},
    success::Bool
)::Path
    wp = Vector{Vector{Float64}}()  # Waypoints
    control_seq = Vector{Vector{Float64}}()  # Control sequence

    if success
        # If goal is reached, extract path to the goal
        wp = [goal]
        current_idx = length(vertices)
    else
        # If goal is not reached, select a random node
        current_idx = rand(1:length(vertices))  # Pick a random node
        println("Goal not reached. Debugging with random node: $current_idx")
        wp = [vertices[current_idx]]
    end

    # Backtrack to start using parents dictionary
    while haskey(parents, current_idx)
        current_idx = parents[current_idx]  # Get the parent index
        push!(wp, vertices[current_idx])
        push!(control_seq, controls[current_idx-1])  # Controls are indexed by edge
    end

    # Reverse the waypoints and controls to maintain correct order
    wp = reverse(wp)
    control_seq = reverse(control_seq)

    # Compute path length
    path_length = sum(norm(wp[i+1] - wp[i]) for i in 1:(length(wp)-1))

    return Path(wp, control_seq, path_length)
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
    best_control = randControl(agent.controlBounds)
    best_score = -Inf
    @distributed for i in 1:25
        u = randControl(agent.controlBounds)
        new_point = propagate(xk, u, dt, agent.EOM, wind)
        score = scoreControl(new_point, xkp1)
        if score > best_score
            best_score = score
            best_control = u
        end
    end
    return best_control
end

function scoreControl(new_point::Vector{Float64}, goal_point::Vector{Float64})::Float64
    x, y, z, v, ψ, γ = new_point[1:6]
    goal_distance_score = -norm(new_point[1:3] - goal_point[1:3])
    ideal_airspeed = 25.0
    airspeed_penalty = abs(v - ideal_airspeed)
    airspeed_score = -airspeed_penalty
    max_glide_slope = π / 4
    attitude_penalty = if abs(γ) > max_glide_slope abs(γ) - max_glide_slope else 0.0 end
    attitude_score = -attitude_penalty

    w_goal, w_airspeed, w_attitude = 1.0, 0.5, 0.5
    return w_goal * goal_distance_score + w_airspeed * airspeed_score + w_attitude * attitude_score
end

function isPointValid(point::Vector{Float64}, lim::Bounds)::Bool
    x, y, z, v, ψ, γ = point[1:6]
    return (lim.xmin <= x <= lim.xmax && lim.ymin <= y <= lim.ymax && lim.zmin <= z <= lim.zmax && 11.0 <= v <= 30.0)
end

function propagate(x::Vector{Float64}, u::Vector{Float64}, dt::Float64, eom::Function, wind::Function)::Vector{Float64}
    return RK4_int_wctrl(dt, x, u, eom, wind)
end

function RK4_int_wctrl(dt, state, u, fdot::Function, wind::Function)::Vector{Float64}
    k1 = dt * fdot(state, u, wind)
    k2 = dt * fdot(state + 0.5 * k1, u, wind)
    k3 = dt * fdot(state + 0.5 * k2, u, wind)
    k4 = dt * fdot(state + k3, u, wind)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end
