using Graphs, NearestNeighbors, Distributed

mutable struct Path
    path::Vector{Vector{Float64}}       # Vector of waypoints, starting from time 0, length t*dt
    control::Vector{Vector{Float64}}    # Vector of controls, starting from time dt, length t*dt-1
    length::Float64                     # Path length in Euclidean space
end
function plan(world::PlanningProblem)::Path
    g = SimpleGraph()
    vertices = Vector{Vector{Float64}}()  # State space points
    parents = Dict{Int,Union{Int,Nothing}}()
    parents[1] = nothing  # Root has no parent
    # initialize
    init_point = world.env.init
    add_vertex!(g)
    push!(vertices, init_point)
    # env knowns
    goal = world.env.goal
    lim = world.env.lim
    max_iter = 1500
    min_dist_thresh = 75.0
    max_outgoing_connections = 10  # limit connections per node
    kdtree_update_frequency = 50  # Update KDTree every 50 iterations
    kdtree = KDTree(hcat([init_point[1:3]]...))

    min_dt, max_dt = 3.0, 12.0

    for i in 1:max_iter
        if i % kdtree_update_frequency == 0
            kdtree = KDTree(hcat([v[1:3] for v in vertices]...))
        end

        # adaptive goal sampling
        goal_sampling_prob = 0.1 + 0.9 * (length(vertices) / max_iter)

        # biased random sample
        rand_point = if rand() < goal_sampling_prob
            goal[1:3]
        else
            biasedSample(lim, kdtree, vertices, goal[1:3])
        end

        nearest_idx = knn(kdtree, rand_point, 1)[1][1]
        nearest_x = vertices[nearest_idx]

        # dynamic timestep based on proximity to the goal
        nearest_dist_to_goal = norm(nearest_x[1:3] - goal[1:3])
        dt = if nearest_dist_to_goal > 300
            max_dt  # Use max timestep for far points
        else
            min_dt + (max_dt - min_dt) * (1 - exp(-0.1 * nearest_dist_to_goal))
        end

        # get best control from set
        ustar = getBestControl(nearest_x, rand_point, dt, world.agents[1], world.env.windField, goal)

        # Propagate dynamics
        new_x = propagate(nearest_x, ustar, dt, world.agents[1].EOM, world.env.windField)

        newnode_dist = norm(new_x[1:3] - nearest_x[1:3])
        # add the new vertex if valid and connection limit is not exceeded
        if isPointValid(new_x, lim) && newnode_dist > min_dist_thresh && degree(g, nearest_idx) < max_outgoing_connections
            add_vertex!(g)
            push!(vertices, new_x)
            add_edge!(g, nearest_idx, length(vertices))
            parents[length(vertices)] = nearest_idx

            #  connect to the goal if close enough
            if norm(new_x[1:3] - goal[1:3]) <= 100
                add_vertex!(g)
                push!(vertices, goal)
                add_edge!(g, length(vertices) - 1, length(vertices))
                parents[length(vertices)]=nearest_idx
                println("Connected to goal!")
                plot_rrt_tree(vertices,parents,lim,goal)
                return extract_path(g, vertices, [[0.0]], parents, init_point, goal, true)
            end
        end
    end
    plot_rrt_tree(vertices,parents,lim,goal)
    return extract_path(g, vertices, [[0.0]], parents, init_point, goal, false)
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
    #control_seq = Vector{Vector{Float64}}()  # Control sequence
    if success
        # If goal is reached, extract path to the goal
        wp = [goal]
        current_idx = length(vertices)
    else
        error("Goal not reached")
    end

    # Backtrack to start using parents dictionary
    while parents[current_idx] !== nothing
        current_idx = parents[current_idx]  # get parent index
        push!(wp, vertices[current_idx])
        #push!(control_seq, controls[current_idx])  
    end

    # Ensure the root (start) is included in the path
    push!(wp, start)

    # Reverse the waypoints and controls to maintain correct order
    wp = reverse(wp)
    #control_seq = reverse(control_seq)

    # Compute path length
    #path_length = sum(norm(wp[i+1] - wp[i]) for i in 1:(length(wp)-1))

    return Path(wp, [[0.0]], 0.0)
end
# sampling functions 
function samplePoint(l::Bounds)::Vector{Float64}
    x = rfinbounds(l.xmin, l.xmax)
    y = rfinbounds(l.ymin, l.ymax)
    z = rfinbounds(l.zmin, l.zmax)
    return [x, y, z]
end
function biasedSample(lim::Bounds, kdtree::KDTree, vertices::Vector{Vector{Float64}}, goal::Vector{Float64})::Vector{Float64}
    for _ in 1:10  # limit max recursion
        rand_point = samplePoint(lim)
        nearest_idx = knn(kdtree, rand_point, 1)[1][1]
        nearest_x = vertices[nearest_idx]
        distance_to_nearest = norm(rand_point - nearest_x[1:3])
        distance_to_goal = norm(rand_point - goal)
        if distance_to_nearest > 50 || rand() < exp(-distance_to_goal / 200.0)
            return rand_point
        end
    end
    return goal  # failsafe
end

# control functions (most of which are self explanatory)
function randControl(ubounds::Vector{Vector{Float64}})::Vector{Float64}
    rf = zeros(Float64, length(ubounds))
    for i in 1:length(ubounds)
        rf[i] = rfinbounds(ubounds[i][1], ubounds[i][2])
    end
    return rf
end
function getBestControl(xk::Vector{Float64}, xkp1::Vector{Float64}, dt::Float64, agent::Agent, wind::Function, goal::Vector{Float64})::Vector{Float64}
    #xk current state, xkp1 proposed sampled state, 
    best_control = randControl(agent.controlBounds)
    best_score = -Inf
    @distributed for i in 1:36  # multithread this (wow)
        u = randControl(agent.controlBounds)
        new_point = propagate(xk, u, dt, agent.EOM, wind)
        score = scoreControl(xk, new_point, xkp1,wind)
        if score > best_score
            best_score = score
            best_control = u
        end
    end
    return best_control
end
function scoreControl(
    old_point::Vector{Float64},
    new_point::Vector{Float64},
    rand_point::Vector{Float64},
    wind::Function
)::Float64
    xkp1, ykp1, zkp1, vkp1, ψkp1, γkp1 = new_point[1:6]
    xk, yk, zk, vk, ψk, γk = old_point[1:6]

    goal_distance_score = -norm(new_point[1:3] - new_point[1:3])  # Lower distance to goal is better

    energy_score = (((vkp1^2) / 2) + zkp1 * 9.81) - (((vk^2) / 2) + zk * 9.81)  # Energy change

    wind_vector = wind(new_point[1:3])
    wind_alignment_score = wind_vector[3]  # Prioritize updraft regions
    # weights
    w_new = 1.0          #  for progress toward goal
    w_energy = 10.0        #  for energy efficiency
    w_wind = 10.0          # Weight for updraft location (also energy, but different just cuz)

    total_score = (
        w_new * goal_distance_score +
        w_energy * energy_score #+
        #w_wind * wind_alignment_score
    )

    return total_score
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