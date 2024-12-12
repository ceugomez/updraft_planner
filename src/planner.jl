using Graphs
using NearestNeighbors
using Distributed

mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # Vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # Path length in Euclidean space
end

function plan(world::PlanningProblem)::Path
    Random.seed!("myproject")
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # Stores state space points corresponding to graph vertices
    controls = Vector{Vector{Float64}}()  # Stores control inputs corresponding to edges
    parents = Dict{Int, Union{Int, Nothing}}()
    parents[1] = nothing  # Root has no parent

    init_point = world.env.init
    add_vertex!(g)
    push!(vertices, init_point)
    
    goal = world.env.goal
    lim = world.env.lim
    max_iter = 10_000
    min_dist_thresh = 75.0
    kdtree_update_frequency = 100  # Update KDTree every 100 iterations
    kdtree = KDTree(hcat([init_point[1:3]]...))  # Initial KDTree

    # Timestep constraints
    dt_min = 1.0
    dt_max = 10.0

    for i in 1:max_iter
        # Update KDTree periodically
        if i % kdtree_update_frequency == 0
            kdtree = KDTree(hcat([v[1:3] for v in vertices]...))
        end

        # Sample a random point
        goal_sampling_prob = min(0.5, length(vertices) / max_iter)
        rand_point = if rand() < goal_sampling_prob
            goal[1:3]
        else
            biasedSample(lim, kdtree, vertices, goal[1:3])
        end

        # Find the nearest vertex
        nearest_idx = knn(kdtree, rand_point, 1)[1][1]
        nearest_x = vertices[nearest_idx]

        # Calculate the dynamic timestep
        dist_to_goal = norm(nearest_x[1:3] - goal[1:3])
        dt = clamp(dist_to_goal / 100.0, dt_min, dt_max)

        # Find the best control
        ustar = getBestControl(nearest_x, rand_point, dt, world.agents[1], world.env.windField, world.env.goal)

        # Propagate dynamics
        new_x = propagate(nearest_x, ustar, dt, world.agents[1].EOM, world.env.windField)

        # Distance check from the nearest node
        newnode_dist = norm(new_x[1:3] - nearest_x[1:3])

        # Add the new vertex if valid
        if isPointValid(new_x, lim) && newnode_dist > min_dist_thresh
            add_vertex!(g)
            push!(vertices, new_x)
            push!(controls, ustar)
            add_edge!(g, nearest_idx, length(vertices))
            parents[length(vertices)] = nearest_idx

            # Check if the goal is reached
            if norm(new_x[1:3] - goal[1:3]) <= 100
                println("Reached goal!")
                plot_rrt_tree(vertices, parents, lim, goal)
                return extract_path(g, vertices, controls, parents, init_point, goal, true)
            end
        end
    end

    println("RRT failed to find a path within the maximum number of iterations")
    plot_rrt_tree(vertices, parents, lim, goal[1:3])
    return extract_path(g, vertices, controls, parents, init_point, goal, false)
end


# backtracking path
function extract_path(
    g::SimpleGraph,
    vertices::Vector{Vector{Float64}},
    controls::Vector{Vector{Float64}},
    parents::Dict{Int,Union{Int,Nothing}},
    start::Vector{Float64},
    goal::Vector{Float64},
    success::Bool,
)::Path
    wp = Vector{Vector{Float64}}()  # Waypoints
    control_seq = Vector{Vector{Float64}}()  # Control sequence
    if success
        # If goal is reached, extract path to the goal
        wp = [goal]
        current_idx = length(vertices)
    else
        # If goal is not reached, select a random node
        current_idx = length(vertices) # Pick a random node
        println("Goal not reached. Debugging with random node: $current_idx")
        wp = [vertices[current_idx]]
    end

    # Backtrack to start using parents dictionary
    while parents[current_idx] !== nothing
        current_idx = parents[current_idx]  # Get the parent index
        push!(wp, vertices[current_idx])
        push!(control_seq, controls[current_idx])  # Controls are indexed by edge
    end

    # Ensure the root (start) is included in the path
    push!(wp, start)

    # Reverse the waypoints and controls to maintain correct order
    wp = reverse(wp)
    control_seq = reverse(control_seq)

    # Compute path length
    path_length = sum(norm(wp[i+1] - wp[i]) for i in 1:(length(wp)-1))

    return Path(wp, control_seq, path_length)
end
# sampling functions 
function samplePoint(l::Bounds)::Vector{Float64}
    x = rfinbounds(l.xmin, l.xmax)
    y = rfinbounds(l.ymin, l.ymax)
    z = rfinbounds(l.zmin, l.zmax)
    return [x, y, z]
end
function biasedSample(
    lim::Bounds,
    kdtree::KDTree,
    vertices::Vector{Vector{Float64}},
    goal::Vector{Float64}
)::Vector{Float64}
    # Generate a random sample point
    rand_point = samplePoint(lim)

    # Find the nearest vertex and compute distances
    nearest_idx = knn(kdtree, rand_point, 1)[1][1]
    nearest_x = vertices[nearest_idx]
    distance_to_nearest = norm(rand_point - nearest_x[1:3])
    distance_to_goal = norm(rand_point - goal)

    # Threshold for considering unexplored regions
    exploration_threshold = 150.0

    # Bias sampling
    if distance_to_nearest > exploration_threshold
        # Accept points in unexplored regions
        return rand_point
    else
        # Accept points near existing nodes with probability inversely proportional to distance from goal
        acceptance_probability = exp(-distance_to_goal / 200.0)  # Tune the scaling factor (200.0) as needed
        if rand() < acceptance_probability
            return rand_point  # Accept the point
        else
            return biasedSample(lim, kdtree, vertices, goal)  # Resample
        end
    end
end
# control functions (most of which are self explanatory)
    function randControl(ubounds::Vector{Vector{Float64}})::Vector{Float64}
        rf = zeros(Float64, length(ubounds))
        for i in 1:length(ubounds)
            rf[i] = rfinbounds(ubounds[i][1], ubounds[i][2])
        end
        return rf
    end
    function getBestControl(xk::Vector{Float64}, xkp1::Vector{Float64}, dt::Float64, agent::Agent, wind::Function,goal::Vector{Float64})::Vector{Float64}
        #xk current state, xkp1 proposed sampled state, 
        best_control = randControl(agent.controlBounds)
        best_score = -Inf
        @distributed for i in 1:150  # multithread this (wow)
            u = randControl(agent.controlBounds)
            new_point = propagate(xk, u, dt, agent.EOM, wind)
            score = scoreControl(xk,new_point,xkp1,goal)
            if score > best_score
                best_score = score
                best_control = u
            end
        end
        return best_control
    end
    function scoreControl(old_point::Vector{Float64}, new_point::Vector{Float64}, rand_point::Vector{Float64}, goal::Vector{Float64})::Float64
        xkp1, ykp1, zkp1, vkp1, ψkp1, γkp1 = new_point[1:6]
        xk, yk, zk, vk, ψk, γk = old_point[1:6]
        # invalidate controls that put us below goal 
        if (xkp1 < goal[3])
            return -Inf
        end
        # close twds goal
        goal_distance_score = -norm(new_point[1:2] - goal[1:2]) # just care about lateral distance to goal
        # gain energy
        energy_score = (((vkp1^2) / 2) + zkp1 * 9.81) - (((vk^2) / 2) + zkp1 * 9.81) # do we gain or lose energy with this control?
        w_goal, w_energy = 1.0, 0.5
        return w_goal * goal_distance_score + w_energy * energy_score
    end
# collision checker (validity, more than anything)
    function isPointValid(point::Vector{Float64}, lim::Bounds)::Bool
        x, y, z, v, ψ, γ = point[1:6]
        if v < 15 || v > 45 || x < lim.xmin || x > lim.xmax || 
        y < lim.ymin || y > lim.ymax || z < lim.zmin || z > lim.zmax
            return false
        end
        return true
    end
# dynamics and integrators
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
# roi selection algorithm
    function select_best_roi(
        current_state::Vector{Float64},
        rois::Vector{Vector{Float64}},
        lim::Bounds,
        agent::Agent,
        wind::Function
    )::Union{Int,Nothing}
        best_cost = Inf
        best_idx = nothing

        for (i, roi) in enumerate(rois)
            energy_cost, feasible = estimate_energy_cost(current_state, roi, agent, wind, lim)
            if feasible && energy_cost < best_cost
                best_cost = energy_cost
                best_idx = i
            end
        end

        return best_idx
    end

# Function to estimate energy cost to an ROI
function estimate_energy_cost(
    current_state::Vector{Float64},
    roi::Vector{Float64},
    agent::Agent,
    wind::Function,
    lim::Bounds
)::Tuple{Float64,Bool}
    num_samples = 10
    energy_cost = 0.0
    feasible = true

    # Sample points along the straight-line path to the ROI
    path_points = [current_state[1:3] .+ (roi - current_state[1:3]) .* (t / num_samples) for t in 0:num_samples]
    for point in path_points
        wind_effect = wind(point)
        cost_increment = compute_energy_increment(point, wind_effect, agent)
        energy_cost += cost_increment

        # Check if the point is within bounds
        if !isPointValid(vcat(point, current_state[4:6]), lim)
            feasible = false
            break
        end
    end
    return (energy_cost, feasible)
end

# Placeholder function for energy increment calculation
function compute_energy_increment(pos::Vector{Float64}, wind_effect::Vector{Float64}, agent::Agent)::Float64
    return norm(wind_effect)  # Replace with more accurate energy computation as needed
end