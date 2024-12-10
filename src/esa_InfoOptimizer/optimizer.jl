using Optimization, PyPlot, BenchmarkTools, LinearAlgebra, Random, Distributions
# structure to store optimization results
struct optimresult
    result::Vector{Float64}
    time::Float64
    iterations::Int64
    fevals::Int64
    fnval::Float64
end
# global structure to store problem statement
struct OptimProblemStatement
    agents::Int64
    dims::Matrix{Int64}
    OptMethod::String
    pos::Vector{Vector{Float64}}
	R::Float64
end
# convenience function - print problem statement 
    function printStartupScript(params::OptimProblemStatement)
        println("Running Targeted Observation Optimization -------------");
        println("cgf cego6160@colorado.edu 11.18.24")
        println("Initializing...")
        println(string(params.agents) * " Agents")
        println("Dimensions: " * string(params.dims[1, 2] - params.dims[1, 1]) * " x " * string(params.dims[2, 2] - params.dims[2, 1]) * " x " * string(params.dims[3, 2] - params.dims[3, 1]))
        println("Optimization: " * params.OptMethod)
        println("Running...")
    end
# Information Field (hartmann)
    function infoField(xx::Vector{Float64})
        # hartmann 3d, scaled
        # Copyright 2013. Derek Bingham, Simon Fraser University (thanks!!)
        x, y, z = xx[1], xx[2], xx[3]
        xmin, xmax = prob.dims[1, :]
        ymin, ymax = prob.dims[2, :]
        zmin, zmax = prob.dims[3, :]
        x_normalized = (x - xmin) / (xmax - xmin)
        y_normalized = (y - ymin) / (ymax - ymin)
        z_normalized = (z - zmin) / (zmax - zmin)
        xvec_normalized = [x_normalized, y_normalized, z_normalized]
        alpha = [1.0, 1.2, 3.0, 3.2]
        A = [[3.0, 10.0, 30.0] [0.1, 10.0, 35.0] [ 3.0, 10, 30] [0.1, 10, 35]];
        P = [[3689, 1170, 2673] [4699, 4387, 7470] [1091, 8732, 5547] [381, 5743, 8828]]
        P = P.*10^(-4)
        outer = 0;
        for ii in 1:4
            inner = 0;
            for jj in 1:3
                xj = xx[jj].*xvec_normalized[jj]/6+0.3
                Aij = A[jj, ii]
                Pij = P[jj, ii]
                inner = inner + Aij*(xj-Pij)^2
            end
            new = alpha[ii] * exp(-inner)
            outer = outer + new;
        end
        return 50*outer;
    end    
# Random value functions
    # get random value in bounds 
        function rfinbounds(lower_bound, upper_bound)::Float64
            return lower_bound + rand() * (upper_bound - lower_bound)
        end
    # get random point within domain
        function randInDomain()
            # Extract bounds from prob.dims
            xmin, xmax = prob.dims[1, :]
            ymin, ymax = prob.dims[2, :]
            zmin, zmax = prob.dims[3, :]

            # Generate random values within each dimension
            x_rand = rfinbounds(xmin, xmax)
            y_rand = rfinbounds(ymin, ymax)
            z_rand = rfinbounds(zmin, zmax)

            # Return the random point as a tuple or array
            return [x_rand, y_rand, z_rand]
        end
# reward with constraints
    function reward(X)
        # extract position
        x, y, z = X[1], X[2], X[3]

        # get value before modification
        value = infoField(X)

        # extract bounds on domain
        xmin, xmax = prob.dims[1, 1], prob.dims[1, 2]
        ymin, ymax = prob.dims[2, 1], prob.dims[2, 2]
        zmin, zmax = prob.dims[3, 1], prob.dims[3, 2]

        # interior point methods
        function log_barrier_inequality(v, threshold, direction)
            diff = v - threshold
            if direction == :greater  
                if diff > 0
                    return -log(diff)
                else
                    return Inf  
                end
            elseif direction == :less  
                diff = threshold - v
                if diff > 0
                    return -log(diff)
                else
                    return Inf 
                end
            end
        end

        # Compute penalties for boundary constraints
        bounds_penalty = 0.0
        bounds_penalty += log_barrier_inequality(x, xmin, :greater)
        bounds_penalty += log_barrier_inequality(x, xmax, :less)
        bounds_penalty += log_barrier_inequality(y, ymin, :greater)
        bounds_penalty += log_barrier_inequality(y, ymax, :less)
        bounds_penalty += log_barrier_inequality(z, zmin, :greater)
        bounds_penalty += log_barrier_inequality(z, zmax, :less)

        # Repulsion penalty to maintain minimum distance of `prob.R` units
        repulsion_penalty = 0.0
        if !isempty(prob.pos)
            for other_agent in prob.pos
                dist = norm(X - other_agent)
                # Add penalty for violating minimum distance
                repulsion_penalty += log_barrier_inequality(dist, prob.R, :greater)
            end
        end

        # Combine the base value, boundary penalty, and repulsion penalty
        J = -value + bounds_penalty + repulsion_penalty

        return J
    end
# Optimizer: cross-entropy (custom)
    function crossentropy(f, x0, ϵ)
        μ = x0                              # Initial mean guess
        xmin, xmax = prob.dims[1, 1], prob.dims[1, 2]
        ymin, ymax = prob.dims[2, 1], prob.dims[2, 2]
        zmin, zmax = prob.dims[3, 1], prob.dims[3, 2]
        σ2 = Diagonal((2*[xmax-xmin, ymax-ymin, zmax-zmin]))           # Initial covariance guess
        i = 1                               # Iterator
        elites = 100                         # Number of elite samples
        genpop = 1000                        # Number of general samples
        maxi = 1000000                        # Iterator upper bound

        # Use BenchmarkTools for precise timing
        timer = @elapsed begin
            # Optimization loop
            while (i < maxi && norm(σ2) > ϵ)
                σ2 += Diagonal(fill(1e-12, length(x0))) # ensure positive definite
                dist = MvNormal(μ, σ2)       # Build distribution
                X = rand(dist, genpop)      # Get samples from distribution
                samples = Float64[]
                for j in 1:genpop
                    push!(samples, f(X[:, j]))
                end
                p = sortperm(samples)       # Sort samples by function value
                elite = X[:, p[1:elites]]   # Select elite samples

                # Update mean and covariance
                μ = mean(elite, dims=2)[:]
                σ2 = cov(elite')
                i += 1
            end
        end

        # Return optimization result
        return optimresult(μ, timer, i, genpop * i, f(μ))
    end
# Functions to visualize the information field
    function visualizeField(dims::Matrix{Int64}, flags::Vector{Bool})
        # Generate a mesh grid (replacement for Python's MeshGrid)
        function generateMeshGrid(xrange::AbstractVector, yrange::AbstractVector)
            X = repeat(xrange', length(yrange), 1)
            Y = repeat(yrange, 1, length(xrange))
            return X, Y
        end

        # Generate the grid for evaluation
        xarray = collect(dims[1, 1]:0.1:dims[1, 2])
        yarray = collect(dims[2, 1]:0.1:dims[2, 2])
        X, Y = generateMeshGrid(xarray, yarray)

        # Helper function to save plots with consistent settings
        function savePlot(fig, filename)
            savefig(fig, filename, dpi=300, bbox_inches="tight")
        end

        # Compute the information field values for the grid
        Z = [infoField([x, y, prob.pos[1][3]]) for (x, y) in zip(X[:], Y[:])]

        for i in 1:length(flags)
            if flags[i]
                plt = figure()

                if i == 1
                    # Contour Plot
                    println("Plotting Contours...")
                    ax = plt.add_subplot(111)
                    contourf(X, Y, reshape(Z, size(X)), cmap="viridis",
                            levels=collect(minimum(Z):((maximum(Z) - minimum(Z)) / 45):maximum(Z)))
                    for j in 1:length(prob.pos)
                        ax.scatter(prob.pos[j][1], prob.pos[j][2], color="red", s=50, label="Point of Interest").set_zorder(5)
                    end
                    colorbar(label="Information Field Value")
                    title("Information Field Visualization")
                    xlabel("X")
                    ylabel("Y")
                    savePlot(plt, "./figures/contour.png")

                elseif i == 2
                    # Info Field Surface
                    println("Plotting Info Field Surface...")
                    ax = plt.add_subplot(111, projection="3d")
                    surf = ax.plot_surface(X, Y, reshape(Z, size(X)), cmap="viridis", vmin=minimum(Z) * 2)
                    colorbar(surf)
                    for j in 1:length(prob.pos)
                        ax.scatter(prob.pos[j][1], prob.pos[j][2], infoField(prob.pos[j]),
                                color="red", s=150, label="Point of Interest").set_zorder(5)
                    end
                    xlabel("X-axis")
                    ylabel("Y-axis")
                    ax.set_zlabel("Z-axis")
                    title("Information Value Map | x = $(round(prob.pos[1][3], digits=2))")
                    savePlot(plt, "./figures/Info_surface.png")

                elseif i == 3
                    # Reward Function Surface
                    println("Plotting Reward Field Surface...")
                    Z_reward = [reward([x, y, 1]) for (x, y) in zip(X[:], Y[:])]
                    ax = plt.add_subplot(111, projection="3d")
                    surf = ax.plot_surface(X, Y, reshape(Z_reward, size(X)), cmap="viridis", vmin=minimum(Z_reward) * 2)
                    colorbar(surf)
                    for j in 1:length(prob.pos)
                        ax.scatter(prob.pos[j][1], prob.pos[j][2], infoField(prob.pos[j]),
                                color="red", s=150, label="Point of Interest").set_zorder(5)
                    end
                    xlim(prob.dims[1, 1], prob.dims[1, 2])
                    ylim(prob.dims[2, 1], prob.dims[2, 2])
                    zlim(minimum(Z_reward) - 5, maximum(Z_reward) + 5)
                    xlabel("X-axis")
                    ylabel("Y-axis")
                    ax.set_zlabel("Z-axis")
                    title("Reward Function Map | x = $(round(prob.pos[1][3], digits=2))")
                    savePlot(plt, "./figures/Reward_surface.png")
                end
            end
        end

        show()
        return nothing
    end
    function plotInfoFieldContours(resolution::Float64)
        """
        Generate (x, y) and (y, z) contour plots of the information field using the global `prob` structure.

        Parameters:
            resolution::Float64 - Spacing between grid points (smaller values yield finer resolution).
        """
        global prob

        # Extract domain bounds from `prob.dims`
        xmin, xmax = prob.dims[1, :]
        ymin, ymax = prob.dims[2, :]
        zmin, zmax = prob.dims[3, :]

        # Generate grid points
        x_vals = collect(xmin:resolution:xmax)
        y_vals = collect(ymin:resolution:ymax)
        z_vals = collect(zmin:resolution:zmax)

        # Compute (x, y) contour
        X, Y = repeat(x_vals', length(y_vals), 1), repeat(y_vals, 1, length(x_vals))
        Z_xy = [infoField([x, y, 1.0]) for (x, y) in zip(X[:], Y[:])]

        # Plot (x, y) contour
        plt.figure()
        contourf(X, Y, reshape(Z_xy, size(X)), cmap="magma", levels=50)
        colorbar(label="Info Field Value")
        title("Information Field Contour: (x, y)")
        xlabel("x")
        ylabel("y")
        # Overlay agent positions (x, y)
        for (i,agent) in enumerate(prob.pos)
            scatter(agent[1], agent[2], color="red", label="Agent $i", s=50)
        end
        legend()
        savefig("./figures/xy_contour.png", dpi=300, bbox_inches="tight")
        println("Saved (x, y) contour plot to './figures/xy_contour.png'.")

        # Compute (y, z) contour
        Y, Z = repeat(y_vals', length(z_vals), 1), repeat(z_vals, 1, length(y_vals))
        Z_yz = [infoField([0.0, y, z]) for (y, z) in zip(Y[:], Z[:])]

        # Plot (y, z) contour
        plt.figure()
        contourf(Y, Z, reshape(Z_yz, size(Y)), cmap="magma", levels=50)
        colorbar(label="Info Field Value")
        title("Information Field Contour: (y, z)")
        xlabel("y")
        ylabel("z")
        # Overlay agent positions (y, z)
        for (i,agent) in enumerate(prob.pos)
            scatter(agent[2], agent[3], color="red", label="Agent $i", s=50)
        end
        legend()
        savefig("./figures/yz_contour.png", dpi=300, bbox_inches="tight")
        println("Saved (y, z) contour plot to './figures/yz_contour.png'.")
        show()
        return nothing
    end
    function plotInfoFieldLineContours(resolution::Float64)
        """
        Generate (x, y), (y, z), and (x, z) contour plots of the information field using the global `prob` structure,
        and overlay agent positions on the contours.

        Parameters:
            resolution::Float64 - Spacing between grid points (smaller values yield finer resolution).
        """
        global prob

        # Extract domain bounds from `prob.dims`
        xmin, xmax = prob.dims[1, :]
        ymin, ymax = prob.dims[2, :]
        zmin, zmax = prob.dims[3, :]

        # Generate grid points
        x_vals = collect(xmin:resolution:xmax)
        y_vals = collect(ymin:resolution:ymax)
        z_vals = collect(zmin:resolution:zmax)

        # Compute (x, y) contour
        X, Y = repeat(x_vals', length(y_vals), 1), repeat(y_vals, 1, length(x_vals))
        Z_xy = [infoField([x, y, 1.0]) for (x, y) in zip(X[:], Y[:])]

        # Plot (x, y) contour
        plt.figure()
        contour_lines = plt.contour(X, Y, reshape(Z_xy, size(X)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Information Field Contour: (x, y)")
        plt.xlabel("x")
        plt.ylabel("y")
        # Overlay agent positions (x, y)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[1], agent[2], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/xy_contour_lines.png", dpi=300, bbox_inches="tight")
        println("Saved (x, y) contour plot to './figures/xy_contour_lines.png'.")

        # Compute (y, z) contour
        Y, Z = repeat(y_vals', length(z_vals), 1), repeat(z_vals, 1, length(y_vals))
        Z_yz = [infoField([0.0, y, z]) for (y, z) in zip(Y[:], Z[:])]

        # Plot (y, z) contour
        plt.figure()
        contour_lines = plt.contour(Y, Z, reshape(Z_yz, size(Y)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Information Field Contour: (y, z)")
        plt.xlabel("y")
        plt.ylabel("z")
        # Overlay agent positions (y, z)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[2], agent[3], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/yz_contour_lines.png", dpi=300, bbox_inches="tight")
        println("Saved (y, z) contour plot to './figures/yz_contour_lines.png'.")

        # Compute (x, z) contour
        X, Z = repeat(x_vals', length(z_vals), 1), repeat(z_vals, 1, length(x_vals))
        Z_xz = [infoField([x, 0.0, z]) for (x, z) in zip(X[:], Z[:])]

        # Plot (x, z) contour
        plt.figure()
        contour_lines = plt.contour(X, Z, reshape(Z_xz, size(X)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Information Field Contour: (x, z)")
        plt.xlabel("x")
        plt.ylabel("z")
        # Overlay agent positions (x, z)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[1], agent[3], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/xz_contour_lines.png", dpi=300, bbox_inches="tight")
        println("Saved (x, z) contour plot to './figures/xz_contour_lines.png'.")

        plt.show()
        return nothing
    end
    function plotRewardFieldLineContours(resolution::Float64)
        """
        Generate (x, y), (y, z), and (x, z) contour plots of the information field using the global `prob` structure,
        and overlay agent positions on the contours.

        Parameters:
            resolution::Float64 - Spacing between grid points (smaller values yield finer resolution).
        """
        global prob

        # Extract domain bounds from `prob.dims`
        xmin, xmax = prob.dims[1, :]
        ymin, ymax = prob.dims[2, :]
        zmin, zmax = prob.dims[3, :]

        # Generate grid points
        x_vals = collect(xmin:resolution:xmax)
        y_vals = collect(ymin:resolution:ymax)
        z_vals = collect(zmin:resolution:zmax)

        # Compute (x, y) contour
        X, Y = repeat(x_vals', length(y_vals), 1), repeat(y_vals, 1, length(x_vals))
        Z_xy = [reward([x, y, 1.0]) for (x, y) in zip(X[:], Y[:])]

        # Plot (x, y) contour
        plt.figure()
        contour_lines = plt.contour(X, Y, reshape(Z_xy, size(X)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Reward Field Contour: (x, y)")
        plt.xlabel("x")
        plt.ylabel("y")
        # Overlay agent positions (x, y)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[1], agent[2], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/xy_contour_lines_reward.png", dpi=300, bbox_inches="tight")
        println("Saved (x, y) contour plot to './figures/xy_contour_lines_reward.png'.")

        # Compute (y, z) contour
        Y, Z = repeat(y_vals', length(z_vals), 1), repeat(z_vals, 1, length(y_vals))
        Z_yz = [infoField([0.0, y, z]) for (y, z) in zip(Y[:], Z[:])]

        # Plot (y, z) contour
        plt.figure()
        contour_lines = plt.contour(Y, Z, reshape(Z_yz, size(Y)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Reward Field Contour: (y, z)")
        plt.xlabel("y")
        plt.ylabel("z")
        # Overlay agent positions (y, z)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[2], agent[3], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/yz_contour_lines_reward.png", dpi=300, bbox_inches="tight")
        println("Saved (y, z) contour plot to './figures/yz_contour_lines_reward.png'.")

        # Compute (x, z) contour
        X, Z = repeat(x_vals', length(z_vals), 1), repeat(z_vals, 1, length(x_vals))
        Z_xz = [reward([x, 0.0, z]) for (x, z) in zip(X[:], Z[:])]

        # Plot (x, z) contour
        plt.figure()
        contour_lines = plt.contour(X, Z, reshape(Z_xz, size(X)), cmap="viridis", levels=20)
        plt.clabel(contour_lines, inline=1, fontsize=8)
        plt.title("Reward Field Contour: (x, z)")
        plt.xlabel("x")
        plt.ylabel("z")
        # Overlay agent positions (x, z)
        for (i, agent) in enumerate(prob.pos)
            scatter(agent[1], agent[3], color="red", label="Agent $i", s=50)
        end
        plt.legend()
        savefig("./figures/xz_contour_lines_reward.png", dpi=300, bbox_inches="tight")
        println("Saved (x, z) contour plot to './figures/xz_contour_lines_reward.png'.")

        plt.show()
        return nothing
    end
# Whole-shebang optimizer function - returns vector of optimized ROIs
    function OptimizeEnvironment(problem::OptimProblemStatement, tol::Float64)
        # OptimProblemStatement(3, [-30.0 30.0; -30.0 30.0; 0.0 2.0], "Cross-Entropy", Vector{Vector{Float64}}(),5.0)
        global prob = problem # Number of agents to consider, Environment dimension, opt method
        printStartupScript(prob)
        x0 = randInDomain() #[0.0,0.0,0.0]#randInDomain()               # Starting guess
        #tol = 1e-5
        # Iteratively solve position problem 
        for i in 1:prob.agents
            # optimize
            result = crossentropy(reward, x0, tol)
            # output position
            println("Agent Position: " * string(trunc(result.result[1], digits=3, base=10)) * ", " * string(trunc(result.result[2], digits=3, base=10)) * ", " * string(trunc(result.result[3], digits=3, base=10)))
            push!(prob.pos, result.result);
            println("Time: $(result.time)");
            println("Fn Evals: $(result.fevals)");
            temp_result = infoField(prob.pos[i])
            println("Fn Value: $(temp_result)");
        end
        out = prob.pos
        prob = Nothing
        return out
    end
function Plot()
    plotInfoFieldContours(0.1)  
    plotInfoFieldLineContours(0.1)
    plotRewardFieldLineContours(0.1)
end
# Main script
# begin
#     prob = OptimProblemStatement(3, [-30.0 30.0; -30.0 30.0; 0.0 2.0], "Cross-Entropy", Vector{Vector{Float64}}(),5.0)
#     OptimizeEnvironment(prob, 1e-5); 
# end
