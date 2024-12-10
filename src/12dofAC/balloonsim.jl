include("generic.jl");
using Parameters, Atmosphere

@with_kw mutable struct balloonParams
    # known from swenson et al; 
    # payload mass: 91g 
    # 125 L volume, ~ 0.125 m^3
    gravity::Float64 = 9.81; 
    ρa::Float64 = 1.225;            # kg⋅m³         (dynamically modified)
    ρh::Float64 = 0.166;            # kg⋅m³         (helium)
    Cd::Float64 = 0.5;              # drag coeff, dynamically modified
    r::Float64 = 0.310175245;       # m            
    vol::Float64 = 0.125;           # m³,          (swenson et.al.)
    m::Float64 = 0.092;             # kg           (swenson et.al)
    AGL::Float64 = 1611.7824;       # m, AGL->MSL conversion  (KBDU)
end
function balloonstatederiv(x::Vector{Float64})
    # Buoyancy force
        rho_z, mu, _ = atmospherefit(x[3]+prm.AGL)
        prm.ρa = rho_z
        fb = (prm.ρa - prm.ρh) * prm.gravity * prm.vol                  
        Fb = [0.0, 0.0, fb]  # Only in z direction

    # Drag force
        W = WindField(x[1:3])
        if (isnan(W[1]) || isnan(W[2]) || isnan(W[3]))
            println("error: NAN in wind field");
        end
        Vrel = x[4:6] - W  # Relative velocity
        # speed 
        rel_speed = norm(Vrel)
        prm.Cd = getDragCoeff(rel_speed,prm.ρa,mu);
        area = π * prm.r^2  
        Fd_mag = -prm.Cd * 0.5 * prm.ρa * rel_speed^2 * area
        Fd = Vrel./rel_speed * Fd_mag   # normlized direction*magnitude
    # Gravity force
        Fg = [0.0, 0.0, -prm.m * prm.gravity]
    # Acceleration: a = ∑F / m
        a = (Fb + Fd + Fg)./prm.m
    # State derivative: [velocity; acceleration]
    xdot = vcat(x[4:6], a)
    return xdot
end
function getDragCoeff(vrel::Float64, rho::Float64, mu::Float64)
    # calculate reynolds number
    RE = rho*abs(vrel)*2*prm.r/mu;
    Cd = 0.0;
    # calculate drag coefficient based on spherical drag model
    if (RE≤1)
        Cd = 24/RE
    end
    if (RE > 1 && RE ≤ 400)
        Cd = 24/(RE^0.646)
    end
    if (RE>400 && RE < 3*10^5)
        Cd = 0.5 
    end
    if (RE>3*10^5)
        Cd = 3.66*10^(-4)*(RE^0.4275)
    end
    return Cd
end
function visualize_field_xyz(history::Matrix{Float64})
    # Validate input dimensions
    if size(history, 1) < 3
        error("Input matrix must have at least 3 rows representing E, N, and Up.")
    end

    # Create subplots for the tracks
    plot1 = Plots.plot(history[1, :]/1000, history[2, :]/1000, xlabel="Easting [m]", ylabel="Northing [m]", title="NE Track")
    plot2 = Plots.plot(history[1, :]/1000, history[3, :], xlabel="Easting [m]", ylabel="Up [m]", title="EU Track")
    plot3 = Plots.plot(history[2, :]/1000, history[3, :], xlabel="Northing [m]", ylabel="Up [m]", title="3D Track")

    # Combine the plots into a grid layout
    layout = @layout [a b; c]
    Plots.plot(plot1, plot2, plot3, layout=layout, size=(900, 600))

    # Save the figure
    Plots.savefig("./figures/Track.png")
    show()
    return nothing
end

function visualize_field_t(history::Matrix{Float64}, dt::Float64, tf::Float64)
    # Validate input dimensions
    @assert size(history, 1) == 3 "Input 'history' must have 3 rows (Easting, Northing, Up)."

    # Create time vector in minutes
    t = LinRange(0, size(history, 2) / (60/dt), size(history, 2))

    # Set up the figure
    plt = figure(figsize=(12, 4))
    plt.suptitle("Balloon Time History", fontsize=14)

    # Plot Easting vs. Time
    ax1 = plt.add_subplot(131)
    ax1.plot(t, history[1, :]/1000, label="Easting", color="blue")
    ax1.set_xlabel("Time [min]")
    ax1.set_ylabel("Easting [km]")
    ax1.grid(true)
    ax1.legend()

    # Plot Northing vs. Time
    ax2 = plt.add_subplot(132)
    ax2.plot(t, history[2, :]/1000, label="Northing", color="green")
    ax2.set_xlabel("Time [min]")
    ax2.set_ylabel("Northing [km]")
    ax2.grid(true)
    ax2.legend()

    # Plot Up vs. Time
    ax3 = plt.add_subplot(133)
    ax3.plot(t, history[3, :], label="Up", color="red")
    ax3.set_xlabel("Time [min]")
    ax3.set_ylabel("Up [m]")
    ax3.grid(true)
    ax3.legend()

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for the title
    PyPlot.savefig("./figures/TTrack.png", dpi=300)
    return nothing
end
# runtime loop
begin
    # Initialize parameter structure
    global prm = balloonParams()

    # Simulation parameters
    tstep = 0.01  # Time step in seconds
    tf = 1000.0 #5500.0  # Total simulation time in seconds
    noIdx = Int64(tf / tstep) + 1  # Number of time steps including initial state

    # Initialize history array (6 state variables: x, y, z, vx, vy, vz)
    history = zeros(Float64, 6, noIdx)

    # Set initial state 
    history[:, 1] = [10.0, 25.0, 70.0, 0.0, 0.0, 0.0]  # [x, y, z, vx, vy, vz]

    # Simulation loop
    for i in 1:(noIdx - 1)
        # Get current state
        stateK = history[:, i]

        # Integrate to find next state
        stateKp1 = RK4_int(tstep, stateK, balloonstatederiv)

        # Store next state in history
        history[:, i + 1] = stateKp1
    end

    println("Simulation complete.")

    # Define the domain and resolution
    domain_x = (-1000, 1000)  # x domain in meters
    domain_y = (-1000, 1000)  # y domain in meters
    domain_z = (0,1500)     # z domain in meters
    resolution = 100         # Number of grid points
    height_z = 40.0;

    # Plot the wind direction
    #plotFancy(domain_x,domain_y,domain_z,history,50.0,resolution,50.0)
    #PleasePleasePlease(domain_x,domain_y,domain_z,history,50.0,resolution,50.0)
    plot_wind_direction_VF(domain_x, domain_z, history, resolution, 0.2)
    #plot_wind_top_down(domain_x, domain_y, height_z, resolution, 0.2)
    plot_wind_top_down_with_path(domain_x, domain_y, history, height_z, resolution, 0.075)
    visualize_field_t(history[1:3, :],tstep,tf)  
    #visualize_field_xyz(history[1:3, :])
    println("Plotted")
end
