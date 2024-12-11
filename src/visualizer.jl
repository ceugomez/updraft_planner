function plotPath3d(h::Path, g::PlanningProblem )
    # Extract x, y, z coordinates from the input
    x = [point[1] for point in h.path]
    y = [point[2] for point in h.path]
    z = [point[3] for point in h.path]
    # Create the trace for the 3D path
    fig = figure(figsize=(12, 9))
    #
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("Easting");
        ax.set_ylabel("Northing");
        ax.set_zlabel("Altitude");
    ax.plot(
        x, y, z, 
        label="3D Trajectory", 
        color="Black", 
        linewidth=1.5, 
        zorder=8  
    )
    start = g.env.init          
    fin = g.env.goal
    ax.scatter(start[1], start[2], start[3], label="Start", color="green");
    ax.scatter(fin[1], fin[2], fin[3], label="End", color="red");
    show()

end