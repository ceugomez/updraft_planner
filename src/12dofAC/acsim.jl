# cgf cego6160@colorado.edu
# 12DOF aircraft simulation in Julia

# Define parameters struct
mutable struct Parameters
    g::Float64 = 9.81    # m/s^2
    rho::Float64 = 1.210 # kg/m
    S::Float64 = 0.6282  # m
    b::Float64 = 3.067   # m
    c::Float64 = 0.208   # m 
    AR::Float64 = 3.067^2 / 0.6282  # dimless
    m::Float64 = 5.74               # kg
    W::Float64 = 5.74 * 9.81        # N

    Ix::Float64 = 14.5939 / (3.2804 * 3.2804) * 4106 / 12^2 / 32.2
    Iy::Float64 = 14.5939 / (3.2804 * 3.2804) * 3186 / 12^2 / 32.2
    Iz::Float64 = 14.5939 / (3.2804 * 3.2804) * 7089 / 12^2 / 32.2
    Ixz::Float64 = 14.5939 / (3.2804 * 3.2804) * 323.5 / 12^2 / 32.2
    # aero derivs
    CDmin::Float64 = 0.0240
    CLmin::Float64 = 0.2052
    K::Float64 = 0.0549
    e::Float64 = 1 / (0.0549 * (3.067^2 / 0.6282) * π)
    CD0::Float64 = 0.0240 + 0.0549 * 0.2052^2
    K1::Float64 = -2 * 0.0549 * 0.2052
    CDpa::Float64 = 0.0240 + 0.0549 * 0.2052^2
    # motor stuff
    Sprop::Float64 = 0.0707
    Cprop::Float64 = 1.0
    kmotor::Float64 = 30.0
    # aero force derivatives
    CL0::Float64 = 0.2219
    Cm0::Float64 = 0.0519
    CY0::Float64 = 0.0
    Cl0::Float64 = 0.0
    Cn0::Float64 = 0.0

    CLalpha::Float64 = 6.196683
    Cmalpha::Float64 = -1.634010
    CLq::Float64 = 10.137584
    Cmq::Float64 = -24.376066
    CLalphadot::Float64 = 0.0
    Cmalphadot::Float64 = 0.0
    
    CYbeta::Float64 = -0.367231
    Clbeta::Float64 = -0.080738
    Cnbeta::Float64 = 0.080613
    CYp::Float64 = -0.064992
    Clp::Float64 = -0.686618
    Cnp::Float64 = -0.039384
    Clr::Float64 = 0.119718
    Cnr::Float64 = -0.052324
    CYr::Float64 = 0.213412
    # control derivatives
    CLde::Float64 = 0.006776
    Cmde::Float64 = -0.06
    CYda::Float64 = -0.000754
    Clda::Float64 = -0.02
    Cnda::Float64 = -0.000078
    CYdr::Float64 = 0.003056
    Cldr::Float64 = 0.000157
    Cndr::Float64 = -0.000856
    # environmental parameters

end
# Define aircraft state struct
struct AcState
    pe::Vector{Vector{Float64}}    # Position in Earth coordinates
    eang::Vector{Vector{Float64}}  # Euler angles
    vbe::Vector{Vector{Float64}}   # Velocity in the body frame
    veang::Vector{Vector{Float64}} # Angular velocity in the body frame
end
# define aircraft control struct
struct AcCtrl 
    da::Float64
    de::Float64
    dt::Float64
    dr::Float64
end
# State derivatives as a function of state, wind, control, and aircraft parameters
function ACstatederiv(state::AcState, control, param::Parameters)
    # state = [ (x y z)::pe ( ϕ θ ψ)::eang (u v w)::Vbe (p q r)::Veang]
    xe = state.pe   # inertial position
    ve = state.vbe  # inertial velocity, body frame 
    body_angles = state.eang
    body_rates = state.veang
    wind = windField(inertial_position_earth)
    rotBodyInertial = Rbe(eang)
    rotInertialBody = transpose(rotBodyInertial);
    # do the easy ones first
        # 

    

    Fa = ACforces(state,control,param)
    
    # tbd
    statedot = AcState(xedot) 
    return statedot
end
# get forces on the aircraft from state and control
function ACforces(state::AcState, u::AcCtrl)

    # aero forces

    

    return nothing
end
# Function to compute the rotation matrix from Euler angles
function Rbe(eang::Vector{Float64})
    # Extract Euler angles
    φ, θ, ψ = eang[1], eang[2], eang[3]  # Roll, pitch, yaw

    # Rotation matrix from body to inertial frame (East-North-Up)
    Rbe = [
        cos(ψ)*cos(θ) - sin(ψ)*sin(φ)*sin(θ)  -sin(ψ)*cos(φ)  cos(ψ)*sin(θ) + sin(ψ)*sin(φ)*cos(θ);
        sin(ψ)*cos(θ) + cos(ψ)*sin(φ)*sin(θ)   cos(ψ)*cos(φ)  sin(ψ)*sin(θ) - cos(ψ)*sin(φ)*cos(θ);
        -cos(φ)*sin(θ)                         sin(φ)         cos(φ)*cos(θ)
    ]
    return Rbe
    
end
# State derivative of a balloon in a wind field (to be implemented)
function balloonstatederiv(state::AcState)
    # Placeholder for balloon state derivative computation
    return state
end

# Returns inertial wind vector as a function of inertial position (E, N, U)
function windField(pos::Vector{Float64})
    W_N = 10  # North wind component
    W_E = 0   # East wind component
    W_U = 1   # Updraft component 
    return [W_N, W_E, W_U]
end

# Runge-Kutta 4th-order integrator
function RK4_int(dt::Float64, state::AcState, fdot)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end
begin
    println("Starting Simulation")
    acparam = Parameters();

end
