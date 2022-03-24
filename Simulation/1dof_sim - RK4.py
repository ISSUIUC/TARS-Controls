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
import src.altimeter as altimeter
import src.kalman_filter as kalman

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

s_dt = .012

# Initialize dictionary to store values at all time steps
ref_a_vels = []

dic = {
       "x":[],
       "x_noise":[],
       "vel": [],
       "accel": [],
       "CD": [],
       "Sref": [],
       "time_sim": [],
       "predict_alt": [],
       "predict_update_alt": []
       }

kalman_dic = {
        "alt": [],
        "vel": []
}

# Initial Values
pos_f = constants.x
pos_f_noise = altimeter.alt_noise(constants.x)
vel_f = constants.vx

# Initial Values for Rk4
rk4_0 = np.array([pos_f, vel_f])
rk4_k = rk4_0

# constructing the acceleration vector 
# ---> derivative of position and velocity
def accel(u, rho, Cd_total):
    r1 = u[0]
    v1 = u[1]
    
    F_a = -((rho* (v1**2) * Sref_a * Cd_total) / 2)
    accel_a = F_a/constants.m0 # Acceleration due to Aerodynamic Forces
    accel_f = accel_a - constants.g

    f1 = np.array([v1, accel_f])

    return f1

def rk4_step(state, dt, rho):
    # rk4 iteration 
    y1 = accel(state, rho)
    y2 = accel(state + 0.5*dt*y1, rho)
    y3 = accel(state + 0.5*dt*y2, rho)
    y4 = accel(state + dt*y3, rho)
    rk4_kp1 = state + dt*(y1 + 2*y2 + 2*y3 + y4)/6
    return rk4_kp1    

def rk4(state, dt):
    
    kalman.initialize(pos_f_noise,vel_f, s_dt)
    
    while (state[1] > 0):
        
        # Define initial flap length at start of control time
        l = 0
        u = np.array([l])
        
        # A-priori 
        # kalman.priori(u)
        
        # grabbing the current states
        pos_f = rk4_k[0]
        pos_f_noise = altimeter.alt_noise(pos_f)
        vel_f = rk4_k[1]
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)
        
        # Total drag coefficient of airframe 
        Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,RASaero,dic["CD"])

        # Approximation - use the area of a circle for reference area
        Sref_a = rocket.sref_approx(constants.D)

        # rk4 iteration 
        rk4_kp1 = rk4_step(state, dt, rho)
        
        # Append Values to the Arrays
        dic["x"].append(float(pos_f))
        dic["x_noise"].append(float(pos_f_noise))
        dic["vel"].append(float(vel_f))
        # dic["accel"].append(float(accel_f?))
        dic["CD"].append(float(Cd_total))
        dic["Sref"].append(float(Sref_a))
        dic["time_sim"].append(float(t))
        
        # A-posteriori update
        kalman.update(pos_f_noise, rk4_kp1[1], Sref_a, rho)
        
        rk4_k = rk4_kp1
                
    return rk4_k
        
        

    
    
    

#* Kalman Filter Initialization
# Initialize states (x), measurement function (H), Covariance [P], White Noise [Q], Measurement Noise Function [R]
kalman.initialize(pos_f_noise,vel_f, s_dt)

# Time setup
start_time = 0
end_time = 30
total_steps = 2500
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
    kalman.priori(u)
    # grabbing the current states
    pos_f = rk4_k[0]
    vel_f = rk4_k[1]
    # Density varies with altitude
    rho = atmosphere.density(pos_f)

    # Total drag coefficient of airframe 
    # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,RASaero,dic["CD"])
    Cd_total = 0.5

    # Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D)

    # rk4 iteration 
    y1 = accel(rk4_k, rho, Cd_total)
    y2 = accel(rk4_k + 0.5*dt*y1, rho, Cd_total)
    y3 = accel(rk4_k + 0.5*dt*y2, rho, Cd_total)
    y4 = accel(rk4_k + dt*y3, rho, Cd_total)
    rk4_kp1 = rk4_k + dt*(y1 + 2*y2 + 2*y3 + y4)/6
    r_kp1 = rk4_kp1[0]; v_kp1 = rk4_kp1[1]

    # calculating acceleration
    # F_a = -((rho* (v_kp1**2) * Sref_a * Cd_total) / 2)
    # accel_a = F_a/constants.m0 # Acceleration due to Aerodynamic Forces
    # accel_f = accel_a - constants.g # Net Acceleration from Aerodynamic Forces + Gravity

    #* End simulation if rocket reached apogee
    if (np.sign(rk4_k[1]) == -1):
        break

    #* Updating Values
    
    rk4_k = rk4_kp1

    # Append Values to the Arrays
    dic["x"].append(float(pos_f))
    dic["x_noise"].append(float(pos_f_noise))
    dic["vel"].append(float(vel_f))
    # dic["accel"].append(float(accel_f?))
    dic["CD"].append(float(Cd_total))
    dic["Sref"].append(float(Sref_a))
    dic["time_sim"].append(float(t))
    

    pos_f_noise = altimeter.alt_noise(pos_f)
    # A-posteriori update
    kalman.update(pos_f_noise, rk4_kp1[1], Sref_a, rho)

    # here we implement another RK4 to find the predicted apogee using the states at this timestep
    # using the previously defined time steps, total_steps can be adjusted based on the amount of computation
    start_time_rk4 = 0
    end_time_rk4 = 30
    total_steps_rk4 = 1000
    time_rk4 = np.linspace(start_time_rk4,end_time_rk4,total_steps_rk4,endpoint=False)
    dt_rk4 = time[1] - time[0]
    # RK4 within the loop starting conditions 
    rk4_0_new = np.array([pos_f, vel_f])
    rk4_k_new = rk4_0_new

    # storing values for this Sim
    x_predicted = []
    v_predicted = []
    #* Starting the RK4 within the loop 
    for t1 in time_rk4:
        # grabbing the current states 
        pos_f_rk4 = rk4_k_new[0]
        vel_f_rk4 = rk4_k_new[1]
        # grabbing current atmospheric properties 
        rho_rk4 = atmosphere.density(pos_f_rk4)
        # Cd_total_rk4 = rasaero.drag_lookup_1dof(pos_f_rk4,vel_f_rk4,RASaero,dic["CD"])
        Cd_total_rk4 = 0.5
        # new rk4 iteration 
        y1_new = accel(rk4_k_new, rho_rk4, Cd_total_rk4)
        y2_new = accel(rk4_k_new + 0.5*dt_rk4*y1_new, rho_rk4, Cd_total_rk4)
        y3_new = accel(rk4_k_new + 0.5*dt_rk4*y2_new, rho_rk4, Cd_total_rk4)
        y4_new = accel(rk4_k_new + dt_rk4*y3_new, rho_rk4, Cd_total_rk4)
        rk4_kp1_new = rk4_k_new + dt_rk4*(y1_new + 2*y2_new + 2*y3_new + y4_new)/6
        r_kp1_new = rk4_kp1_new[0]; v_kp1_new = rk4_kp1_new[1]

        #* End simulation if rocket reached apogee
        if (np.sign(rk4_k_new[1]) == -1):
            break
        # updating the values
        rk4_k_new = rk4_kp1_new
        # parsing the values
        x_predicted.append(pos_f_rk4)
        v_predicted.append(pos_f_rk4)


    
    kalman_dic["alt"].append(kalman.x_k[0][0])
    kalman_dic["vel"].append(kalman.x_k[0][1])
    dic["predict_alt"].append(max(x_predicted))
    
#Print Apogee and total time taken
print("APOGEE (ft):", conversion.m_to_ft(max(dic["x"])))
print("Total Time Taken (s):", t)

#* --------------------------------- Plotting --------------------------------- #
#? Calculating the difference between the predicted altitude and actual altitude 
# difference_pre = []
# difference_post = []
# for num in np.arange(0,len(dic["predict_alt"])):
#     D = dic["predict_alt"][num] - max(dic["x"])
#     d = dic["predict_update_alt"][num] - max(dic["x"])
#     difference_pre.append(D)
#     difference_post.append(d)



# Measurements vs Kalman Filter Graph
plt.plot(dic["time_sim"], dic["x_noise"],label="Noisy Altitude Measurement",color="lightsteelblue",linestyle=":")
plt.plot(dic["time_sim"], dic["x"],label="True Altitude",color="royalblue", linewidth = 3); 
plt.plot(dic["time_sim"], kalman_dic["alt"],label="Estimation",linestyle="--",color="tab:red")
plt.plot(dic["time_sim"], dic["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)
# plt.subplot(1,2,1); plt.plot(dic["time_sim"], dic["predict_alt"],label="Predicted Apogee - Energy Method",linestyle="--", color="tab:green", linewidth = 4.5); plt.legend(fontsize = 10); 
# plt.subplot(1,2,1); plt.plot(dic["time_sim"], dic["predict_update_alt"], label="Corrected Prediction", color="tab:cyan", linestyle="dotted", linewidth = 4.5);
# plt.subplot(1,2,1); 
plt.axhline(y = max(dic["x"]), color = "tab:red", linestyle = "dotted", linewidth = 4.5, label="True Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)

# plt.plot(dic["time_sim"][:-1], difference, label="Difference between Alt_predicted and True", color="tab:blue", linewidth = 3.5, linestyle = "dotted")
# plt.subplot(1,2,2); plt.plot(dic["vel"], difference_pre, label="Pre-correction Error", color="tab:orange", linestyle="dotted", linewidth = 4.5)
# plt.subplot(1,2,2); plt.plot(dic["vel"], difference_post, label="Post-correction Error", color="tab:green", linestyle="dotted", linewidth = 4.5)

plt.xlabel("Time (s)", fontsize = 14)
plt.ylabel("Altitude (m)", fontsize = 14)
plt.legend(fontsize = 14)
plt.show()




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
# plt.plot(time_flight,dic["x"])
# plt.plot(time_flight,dic["x_noise"])
# plt.plot(time_flight,kalman_dic["alt"])
# plt.ylabel("Altitude $(m)$")
# plt.xlabel("Time (s)")
# plt.show()

# # CD Plot
# plt.plot(time_flight,dic["CD"])
# plt.ylabel("CD")
# plt.xlabel("Time (s)")
# plt.show()