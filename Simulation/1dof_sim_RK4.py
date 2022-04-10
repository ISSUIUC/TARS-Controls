from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd
from filterpy.common import Q_continuous_white_noise
import time as timer

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation
import src.RASAero_lookup as rasaero
import src.OpenRocket_lookup as ork
import src.altimeter as altimeter
import src.kalman_filter as kalman
from src.system_propagation import rk4_sim

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
RASaero = pd.read_csv("Simulation/Lookup/RASAero_Intrepid_5800_mk6.csv")
ORK = pd.read_csv("Simulation/OpenRocket Simulations/Intrepid_mk6_April.csv")

April_apogee_time = 32.023  # sec
IREC_apogee_time = 50.041   # sec

April_apogee_goal = 15000  # ft
IREC_apogee_goal = 30000   # ft

April_burnout_alt = 958.75 # m 
IREC_burnout_alt = 1462.25 # m

April_m = 19.0586
IREC_m = 21.1066

April_accel = -31.593 #m/s^2
IREC_accel = -47.825 #m/s^2

constants.m0 = April_m

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I, c_m, m = rocket.I_new(0,0)

# Cd vs Mach Number Polyfit
poly = rasaero.drag_lookup_curve_fit_poly()

# Initial + Desired Values
constants.x = April_burnout_alt
constants.vx = ork.alt_vel_poly_fit(constants.x, ORK, apogee_time=April_apogee_time)
constants.ax = April_accel
init_state = np.array([constants.x, constants.vx])
des_apogee = conversion.ft_to_m(April_apogee_goal) #meters

# Time between sensor readings / KF updates
s_dt = 0.03
# Simulation step-size
dt = 0.006

# Run Sim with and without control
flight_time_nc, kalman_dict_nc, sim_time_nc, sim_dict_nc = rk4_sim(init_state, dt, RASaero, poly, des_apogee, constants.ax)
print("No Control Sim Finished")
flight_time_c, kalman_dict_c, sim_time_c, sim_dict_c = rk4_sim(init_state, dt, RASaero, poly, des_apogee, constants.ax, control=1)

#Print Housekeeping Values
print("APOGEE (No Control) (ft):", conversion.m_to_ft(max(sim_dict_nc["x"])))
print("Flight Time (No Control)(s):", flight_time_nc)
print("Simulator Runtime (No Control)(s): ", sim_time_nc)

print("APOGEE (Control) (ft):", conversion.m_to_ft(max(sim_dict_c["x"])))
print("Flight Time (Control) (s):", flight_time_c)
print("Simulator Runtime (Control) (s): ", sim_time_c)

#* --------------------------------- Plotting --------------------------------- #
#? Calculating the difference between the predicted altitude and actual altitude 
# difference_pre = []
# difference_post = []
# for num in np.arange(0,len(dic["predict_alt"])):
#     D = dic["predict_alt"][num] - max(dic["x"])
#     d = dic["predict_update_alt"][num] - max(dic["x"])
#     difference_pre.append(D)
#     difference_post.append(d)

#* Compare Control vs No Control
plt.subplot(1,2,1)
plt.plot(sim_dict_nc["time_sim"], sim_dict_nc["x"],label="Altitude (No Control)",color="royalblue", linewidth = 3); 
plt.plot(sim_dict_c["time_sim"], sim_dict_c["x"],label="Altitude (Control)",color="green", linewidth = 3); 
plt.plot(sim_dict_nc["time_sim"], sim_dict_nc["x_noise"],label="Noisy Altitude Measurement (No Control)",color="lightsteelblue",linestyle=":")
plt.plot(kalman_dict_nc["time"], kalman_dict_nc["alt"],label="Estimation (No Control)",linestyle="--",color="tab:red")
plt.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
plt.plot(sim_dict_c["time_sim"], sim_dict_c["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)

# plt.subplot(1,2,2)
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["flap_extension"],label="Flap Extension (Control)",color="royalblue", linewidth = 3); 

#* Measurements vs Kalman Filter Graph
# plt.plot(sim_dict_nc["time_sim"], sim_dict_nc["x_noise"],label="Noisy Altitude Measurement",color="lightsteelblue",linestyle=":")
# plt.plot(sim_dict_nc["time_sim"], sim_dict_nc["x"],label="True Altitude",color="royalblue", linewidth = 3); 
# plt.plot(kalman_dict_nc["time"], kalman_dict_nc["alt"],label="Estimation",linestyle="--",color="tab:red")
# plt.axhline(y = max(sim_dict_nc["x"]), color = "tab:red", linestyle = "dotted", linewidth = 4.5, label="True Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
# plt.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)

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