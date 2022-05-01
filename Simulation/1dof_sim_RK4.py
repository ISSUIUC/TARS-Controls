from platform import mac_ver
from turtle import pos
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
import src.propellant_mass as prop

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

# Set to 1 for IREC launch, 0 for April Launch
launch_arg = 1

# RASAero File: Stays the same
# RASaero = pd.read_csv("Simulation/Lookup/RASAero_Intrepid_5800_mk6.csv")
RASaero = pd.read_csv("Simulation/Lookup/RASAero_noAoA_trimmed.csv")

# Importing RasAero Package for Coeffiecient of Drag Lookup
thrust_file = "Simulation/Lookup/AeroTech_M2500T_Trimmed.csv"
constants.apogee_goal = constants.April_apogee_goal
constants.m0 = constants.m0_April
constants.mf = constants.mf_April
constants.thrust_start = constants.thrust_start_April
constants.thrust_end = constants.thrust_end_April
prop_mass_func = prop.find_prop_mass_april

# Change values if simulating IREC launch
if (launch_arg):
    thrust_file = "Simulation/Lookup/Cesaroni_20146N5800-P_Trimmed.csv"
    constants.apogee_goal = constants.IREC_apogee_goal
    constants.m0 = constants.m0_IREC
    constants.mf = constants.mf_IREC
    constants.thrust_start = constants.thrust_start_IREC
    constants.thrust_end = constants.thrust_end_IREC
    prop_mass_func = prop.find_prop_mass_irec
    
thrust_csv = pd.read_csv(thrust_file)

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I, c_m, m = rocket.I_new(0,0)

# Cd vs Mach Number Polyfit
poly_nothrust = rasaero.drag_lookup_curve_fit_poly(0, RASaero)
poly_thrust = rasaero.drag_lookup_curve_fit_poly(1, RASaero)

# Initial + Desired Values
# Position, Velocity, Acceleration are all 0 at launch
init_state = np.array([constants.x, constants.vx, constants.ax])
des_apogee = conversion.ft_to_m(constants.apogee_goal) #meters

# Time between sensor readings / KF updates
s_dt = 0.03
# Simulation step-size
dt = 0.006

# Run Sim with and without control
# flight_time_nc, kalman_dict_nc, sim_time_nc, sim_dict_nc = rk4_sim(init_state, dt, RASaero, poly_nothrust, poly_thrust, des_apogee, thrust_csv, prop_mass_func)
# print("No Control Sim Finished")
flight_time_c, kalman_dict_c, sim_time_c, sim_dict_c = rk4_sim(init_state, dt, RASaero, poly_nothrust, poly_thrust, des_apogee, thrust_csv, prop_mass_func, control=1)

#Print Housekeeping Values
# print("APOGEE (No Control) (ft):", conversion.m_to_ft(max(sim_dict_nc["x"])))
# print("Flight Time (No Control)(s):", flight_time_nc)
# print("Simulator Runtime (No Control)(s): ", sim_time_nc)

print("APOGEE (Control) (ft):", conversion.m_to_ft(max(sim_dict_c["x"])))
print("Flight Time (Control) (s):", flight_time_c)
print("Simulator Runtime (Control) (s): ", sim_time_c)

#* --------------------------------- Plotting --------------------------------- #

#* Compare Control vs No Control
# plt.subplot(1,2,1)
# plt.plot(sim_dict_nc["time_sim"], sim_dict_nc["x"],label="Altitude (No Control)",color="royalblue", linewidth = 3); 
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["x"],label="Altitude (Control)",color="green", linewidth = 3);
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)
# plt.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
# plt.ylabel("Altitude (m)", fontsize = 14)

#* Control Graph (Noisy Measurements, Estimated Altitude, Estimated Apogee)
# plt.subplot(1,2,1)
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["x_noise"],label="Noisy Altitude Measurement",color="lightsteelblue",linestyle=":")
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["x"],label="True Altitude",color="royalblue", linewidth = 3); 
# plt.plot(kalman_dict_c["time"], kalman_dict_c["alt"],label="Estimation",linestyle="--",color="tab:red")
# plt.plot(sim_dict_c["time_sim"], sim_dict_c["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)
# plt.axhline(y = max(sim_dict_c["x"]), color = "tab:red", linestyle = "dotted", linewidth = 4.5, label="True Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
# plt.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
# plt.ylabel("Altitude (m)", fontsize = 14)

# fig,(alt_nc,vel_nc,accel_nc,flap_nc) = plt.subplots(4,1,figsize=(15,10), sharex=True)

# # Altitude Measurements vs Real Altitude vs Kalman Filter Graph (Control)
# alt_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)
# alt_nc.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
# alt_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["x_noise"],label="Noisy Altitude Reading",color="lightsteelblue", linewidth = 3, linestyle=":");
# alt_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["x"],label="Real Altitude",color="royalblue", linewidth = 3);  
# alt_nc.plot(kalman_dict_nc["time"], kalman_dict_nc["alt"],label="Altitude Estimation",linestyle="--",color="tab:red")
# alt_nc.set(ylabel = "Altitude (m)")
# alt_nc.legend()

# # Real Velocity vs Kalman Filter Graph (Control)
# vel_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["vel"],label="Real Velocity",color="royalblue", linewidth = 3);  
# vel_nc.plot(kalman_dict_nc["time"], kalman_dict_nc["vel"],label="Velocity Estimation",linestyle="--",color="tab:red")
# vel_nc.set(ylabel = "Velocity (m/s)")
# vel_nc.legend()

# # Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (Control)
# accel_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["accel_noise"],label="Noisy Accelerometer Reading",color="lightsteelblue", linewidth = 3, linestyle=":");
# accel_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["accel"],label="Real Acceleration",color="royalblue", linewidth = 3);  
# accel_nc.plot(kalman_dict_nc["time"], kalman_dict_nc["accel"],label="Acceleration Estimation",linestyle="--",color="tab:red")
# accel_nc.set(ylabel = "Acceleration (m/s^2)")
# accel_nc.legend()

# flap_nc.plot(sim_dict_nc["time_sim"], sim_dict_nc["flap_extension"],label="Flap Extension (Control)",color="royalblue", linewidth = 3); 
# flap_nc.set(ylabel = "Flap Extension Length (m)")
# flap_nc.legend()


fig,(alt_c,vel_c,accel_c,flap_c) = plt.subplots(4,1,figsize=(15,10), sharex=True)

# Altitude Measurements vs Real Altitude vs Kalman Filter Graph (Control)
alt_c.plot(sim_dict_c["time_sim"], sim_dict_c["predict_alt"], label="Predicted Apogee", linestyle="dashed", color="tab:green", linewidth = 3.5)
alt_c.axhline(y = des_apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Desired Apogee");plt.legend(fontsize = 14); plt.xlabel("Time (s)", fontsize = 14)
alt_c.plot(sim_dict_c["time_sim"], sim_dict_c["x_noise"],label="Noisy Altitude Reading",color="lightsteelblue", linewidth = 3, linestyle=":");
alt_c.plot(sim_dict_c["time_sim"], sim_dict_c["x"],label="Real Altitude",color="royalblue", linewidth = 3);  
alt_c.plot(kalman_dict_c["time"], kalman_dict_c["alt"],label="Altitude Estimation",linestyle="--",color="tab:red")
alt_c.set(ylabel = "Altitude (m)")
alt_c.legend()

# Real Velocity vs Kalman Filter Graph (Control)
vel_c.plot(sim_dict_c["time_sim"], sim_dict_c["vel"],label="Real Velocity",color="royalblue", linewidth = 3);  
vel_c.plot(kalman_dict_c["time"], kalman_dict_c["vel"],label="Velocity Estimation",linestyle="--",color="tab:red")
vel_c.set(ylabel = "Velocity (m/s)")
vel_c.legend()

# Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (Control)
accel_c.plot(sim_dict_c["time_sim"], sim_dict_c["accel_noise"],label="Noisy Accelerometer Reading",color="lightsteelblue", linewidth = 3, linestyle=":");
accel_c.plot(sim_dict_c["time_sim"], sim_dict_c["accel"],label="Real Acceleration",color="royalblue", linewidth = 3);  
accel_c.plot(kalman_dict_c["time"], kalman_dict_c["accel"],label="Acceleration Estimation",linestyle="--",color="tab:red")
accel_c.set(ylabel = "Acceleration (m/s^2)")
accel_c.legend()

flap_c.plot(sim_dict_c["time_sim"], sim_dict_c["flap_extension"],label="Flap Extension (Control)",color="royalblue", linewidth = 3); 
flap_c.set(ylabel = "Flap Extension Length (m)")
flap_c.legend()

fig.tight_layout()
plt.xlabel("Time (s)", fontsize = 14)
plt.show()