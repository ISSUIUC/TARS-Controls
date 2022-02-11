from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd
from filterpy.common import Q_continuous_white_noise

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation
import src.RASAero_lookup as rasaero

#* ---------------------------- Frames we are using --------------------------- #
# Fixed Frame - fixed to the launch rail, not taking into account rotation of the Earth
# X points vertically
# Y is chosen through the right-hand rule
# Z points towards the launch rail

# Body-Fixed Frame - Origin at the center of mass of the rocket
# X points through the nose of the rocket
# Y is chosen through the right-hand rule
# Z points towards the launch lugs on the rocket

# Aerodynamic Frame (Body-Fixed) with Origin at the center-of-mass
# X aligned with aerodynamic velocity
# Y is chosen through the right-hand rule
# Z lies in plane that goes through the centerline of the rocket and the launch lug

#* ---------------------------- Naming Conventions ---------------------------- #
# pos = position(x,y,z)
# or = orientation(yaw,pitch,roll)
# vel = velocity(x,y,z)
# angvel = angular velocity(yaw,pitch,roll)
# accel = acceleration(x,y,z)
# angaccel = angular acceleration (yaw,pitch,roll)

#f = fixed frame
#b = body frame
#a = aerodynamic frame

# pos_f-> position in the fixed frame
# angvel_b-> angular velocity in the body frame
#* ------------------------------ Simulation Code ----------------------------- #

# Importing RasAero Package for Coeffiecient of Drag Lookup
RASaero = pd.read_csv("Simulation/Lookup/RASAero.csv")

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I, c_m, m = rocket.I_new(0,0)

# Initialize dictionary to store values at all time steps
ref_a_vels = []

dic = {
       "x":[],
       "vel": [],
       "accel": [],
       "CD": [],
       "Sref": []
       }

kalman_dic = {
        "alt": [],
        "vel": []
}

# Initial Values
pos_f = constants.x
vel_f = constants.vx

#* Kalman Filter Initialization
# Initialize states (x), measurement function (H)
s_dt = 0.003
x_k = np.array([[pos_f],
                [vel_f]])

# F = np.array([[1. , 0.003],
            #   [0., 0.99999991]])

F = np.array([[1. , 0.003],
              [0., 1]])

H = np.array([float(1),0])
B = np.array([[float(1)],
              [0]])

# Initialize belief in state 
# (Covariance [P], White Noise [Q], Measurement Noise Function [R])
P_k = np.array([[0.01,0],
              [0,1]])
Q = Q_continuous_white_noise(dim=2, dt=s_dt ,spectral_density=1)
R = np.array([1350.0])

# Time setup
start_time = 0
end_time = 30
total_steps = 10000
time = np.linspace(start_time,end_time,total_steps,endpoint=False)
dt = time[1] - time[0]

# Define max and min values for flap actuation
l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

# Define initial flap length at start of control time
l = 0
u = np.array([l])
for t in time:
    
    # A-priori 
    #x(-) = Fx + Bu
    #P(-) = FPF.T + Q
    x_priori = (F @ x_k) + (B @ u)
    P_priori = (F @ P_k @ F.T) + Q
    
    # Density varies with altitude
    rho = atmosphere.density(pos_f)

    # Total drag coefficient of airframe 
    # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,RASaero,dic["CD"])
    Cd_total = 0

    #Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D)

    #Calculate aerodynamic forces on the rocket and the acceleration in the aerodynamic frame
    F_a = -((rho* (vel_f**2) * Sref_a * Cd_total) / 2)
    accel_a = F_a/constants.m0 # Acceleration due to Aerodynamic Forces
    accel_f = accel_a - constants.g # Net Acceleration from Aerodynamic Forces + Gravity

    #* End simulation if rocket reached apogee
    if (np.sign(vel_f) == -1):
        break

    #* Updating Values
    
    # Append Values to the Arrays
    dic["x"].append(float(pos_f))
    dic["vel"].append(float(vel_f))
    dic["accel"].append(float(accel_f))
    dic["CD"].append(float(Cd_total))
    dic["Sref"].append(float(Sref_a))

    # Calculate new velocities and positions using current values
    pos_f = pos_f + (vel_f * dt) + (0.5 * (accel_f * (dt**2)))
    vel_f = vel_f + accel_f*dt
    
    # Received Sensor Measurements
    # A-posteriori update
    # Kalman Gain, posteriori state, Covariance update
    # Update State Guess
    K = P_priori @ H.T * np.reciprocal(H @ P_priori @ H.T + R)
    x_k = x_priori + K @ (np.array([[pos_f],[vel_f]]) - H @ x_priori)
    P_k = (np.eye(2) - K@H) @ P_priori
    
    # F[1][1] = 1 + (Sref_a*rho*0.58*vel_f * s_dt)
    
    kalman_dic["alt"].append(x_k[0][0])
    kalman_dic["vel"].append(x_k[0][1])

    # F[2,2] = 
    
#Print Apogee and total time taken
print("APOGEE (ft):", conversion.m_to_ft(max(dic["x"])))
print("Total Time Taken (s):", t)

#* --------------------------------- Plotting --------------------------------- #
#Calculate the number of steps simulated and create a new linspace
simulated_steps = int(total_steps * ((t+dt - start_time) / (end_time - start_time)))
time_flight = np.linspace(start_time,t,simulated_steps,endpoint=False)


# # Acceleration Plot
# plt.plot(time_flight,dic["accel"])
# plt.ylabel("Acceleration $(m / s^2)$")
# plt.xlabel("Time (s)")
# plt.show()

# # Velocity Plot
# plt.plot(time_flight,dic["vel"])
# plt.ylabel("Vertical Velocity $(m / s)$")
# plt.xlabel("Time (s)")
# plt.show()

# Altitude Plot
plt.plot(time_flight,dic["x"])
plt.plot(time_flight,kalman_dic["alt"])
plt.ylabel("Altitude $(m)$")
plt.xlabel("Time (s)")
plt.show()

# # CD Plot
# plt.plot(time_flight,dic["CD"])
# plt.ylabel("CD")
# plt.xlabel("Time (s)")
# plt.show()