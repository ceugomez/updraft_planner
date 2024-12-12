function plotPath3d(h::Path, g::PlanningProblem)
    # Extract x, y, z coordinates from the input
    x = [point[1] for point in h.path]
    y = [point[2] for point in h.path]
    z = [point[3] for point in h.path]
    # Create the trace for the 3D path
    fig = figure(figsize=(12, 9))
    #
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("X");
        ax.set_ylabel("Y");
        ax.set_zlabel("Z");
    ax.plot(
        x, y, z, 
        label="3D Trajectory", 
        color="Black", 
        linewidth=1.5#, 
        #zorder=8  
    )
    start = g.env.init          
    fin = g.env.goal
    ax.scatter(start[1], start[2], start[3], label="Start", color="green");
    ax.scatter(fin[1], fin[2], fin[3], label="End", color="red");
    ax.scatter(x,y,z, label="Waypoints", color="blue")
    ax.set_xlim([g.env.lim.xmin,g.env.lim.xmax]);
    ax.set_ylim([g.env.lim.ymin,g.env.lim.ymax]);
    ax.set_zlim([g.env.lim.zmin,g.env.lim.zmax])
    
    show()

end
function visualizeWind()
    x = y = z = -100:10:100
    points = [(xi, yi, 50.0) for xi in x, yi in y]

    # Compute wind vectors
    winds = [windField_planner([px, py, pz]) for (px, py, pz) in points]

    # Extract components for quiver plot
    wx = [w[1] for w in winds]
    wy = [w[2] for w in winds]
    wz = [w[3] for w in winds]

    # Quiver plot
    quiver([px for (px, _, _) in points], 
        [py for (_, py, _) in points], 
        quiver=(wx, wy), scale=0.1, color=:blue);
end

function plot_rrt_tree(vertices::Vector{Vector{Float64}}, parents::Dict{Int, Union{Int, Nothing}}, g::Bounds, xg::Vector{Float64})
    fig = figure(figsize=(10, 8))  # Create a new figure
    ax = fig.add_subplot(111, projection="3d")  # 3D subplot

    # Loop through vertices and their parents to plot edges
    for (child_idx, parent_idx) in parents
        if parent_idx !== nothing
            child = vertices[child_idx][1:3]  # Extract x, y, z
            parent = vertices[parent_idx][1:3]  # Extract x, y, z

            # Plot edge between parent and child
            ax.plot(
                [parent[1], child[1]],
                [parent[2], child[2]],
                [parent[3], child[3]],
                color="blue",
                linewidth=0.5
            )
        end
    end

    # Add start and goal points
    start_x, start_y, start_z = vertices[1][1:3]  # Start point
    ax.scatter([start_x], [start_y], [start_z], color="green", s=50, label="Start")  # Start point

    # Goal point
    ax.scatter([xg[1]], [xg[2]], [xg[3]], color="red", s=50, label="Goal")  # Goal point

    # Set labels and title
    ax.set_title("RRT Tree")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_xlim([g.xmin,g.xmax]);
    ax.set_ylim([g.ymin,g.ymax]);
    ax.set_zlim([g.zmin,g.zmax])
    ax.legend()

    # Show the plot
    show()
end
function plot_dubins_glider_path(path::Path, dt::Float64)
    # Extract time vector
    num_states = length(path.path)
    t = [i * dt for i in 0:num_states - 1]

    # Transpose path.path for easier indexing by state
    state_matrix = hcat(path.path...)

    # Extract each state
    x = state_matrix[1, :]
    y = state_matrix[2, :]
    z = state_matrix[3, :]
    v_i = state_matrix[4, :]
    ψ = state_matrix[5, :]
    γ = state_matrix[6, :]

    # Create a 3D plot for trajectory
    fig = figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("Easting")
    ax.set_ylabel("Northing")
    ax.set_zlabel("Altitude")
    ax.plot(x, y, z, label="3D Trajectory", color="Black", linewidth=1.5)
    ax.legend()
    title("3D Trajectory of Dubins Glider")

    # Create individual 2D plots for each state
    fig, axes = subplots(3, 2, figsize=(15, 10))
    axes[1, 1].plot(t, x, label="x (m)", color="blue")
    axes[1, 1].set_title("x vs Time")
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("x (m)")

    axes[1, 2].plot(t, y, label="y (m)", color="orange")
    axes[1, 2].set_title("y vs Time")
    axes[1, 2].set_xlabel("Time (s)")
    axes[1, 2].set_ylabel("y (m)")

    axes[2, 1].plot(t, z, label="z (m)", color="green")
    axes[2, 1].set_title("z vs Time")
    axes[2, 1].set_xlabel("Time (s)")
    axes[2, 1].set_ylabel("z (m)")

    axes[2, 2].plot(t, v_i, label="v_i (m/s)", color="red")
    axes[2, 2].set_title("v_i vs Time")
    axes[2, 2].set_xlabel("Time (s)")
    axes[2, 2].set_ylabel("v_i (m/s)")

    axes[3, 1].plot(t, ψ, label="ψ (rad)", color="purple")
    axes[3, 1].set_title("ψ vs Time")
    axes[3, 1].set_xlabel("Time (s)")
    axes[3, 1].set_ylabel("ψ (rad)")

    axes[3, 2].plot(t, γ, label="γ (rad)", color="brown")
    axes[3, 2].set_title("γ vs Time")
    axes[3, 2].set_xlabel("Time (s)")
    axes[3, 2].set_ylabel("γ (rad)")

    tight_layout()
    show()
end