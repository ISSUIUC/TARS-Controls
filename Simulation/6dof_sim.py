from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation

#TODO: Make dictionary appending into functions
#TODO: Add other values into dictionary
#TODO: Plotting Functions
    # Alpha, Beta, Sref, Yaw, Pitch, Roll
#TODO: Move csv drag function into src library
#TODO: Double check moment of Inertia values

#* Importing RasAero Package
rasaero = pd.read_csv("Lookup/RASAero.csv")
# extracting the columns of interest
mach_num = rasaero.mach.values; aoa = rasaero.alpha_deg.values; cd = rasaero.cd_power_off.values; protub = rasaero.protuberance.values
# narrowing down the columns using mach number range (0.02 - 1.01)
min_index = min(np.where(mach_num == 0.01)[0])
max_index = min(np.where(mach_num == 1.01)[0])
# re-make the lists using this range
mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]

def drag_from_csv(z, velocity_body):
    mach = round(np.linalg.norm(velocity_body) / atmosphere.speed_sound(z),2)
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[1][0]
    vz_b = velocity_body[2][0]
    alpha = abs(round(np.rad2deg(np.arctan2(vz_b,vx_b)), 0))

    mach_index_array = np.where(mach_num == mach)[0]; alpha_index_array = np.where(aoa == alpha)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))

    if len(intersection) == 0:
        return Cd_list[-1]

    index = intersection[0]
    return cd[index]

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

# NumPy arrays to store current states
pos_f = constants.init_pos_f

or_f = constants.init_or_f

vel_f = constants.init_vel_f

angvel_f = constants.init_angvel_f

# Initialize dictionary to store values at all time steps
ref_a_vels = []

dic = {
       "pos_vals":     {"x":[], "y":[], "z":[]},
       "or_vals":      {"Yaw":[], "Pitch":[], "Roll":[]},
       "vel_vals":     {"Vx":[], "Vy":[], "Vz":[]},
       "angvel_vals":  {"Yaw rate":[], "Pitch rate":[], "Roll rate":[]},
       "accel_vals":   {"Ax":[], "Ay":[], "Az":[]},
       "angaccel_vals":{"Yaw accel":[], "Pitch accel":[], "Roll accel":[]}
       }

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
l1 = 0
l2 = 0

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I, c_m, m = rocket.I_new(0,0)
Cd_list = []
for t in time:

    #* Velocity Conversions
    # Rotation from fixed to body frame, Velocity converstion to the body frame
    R_fb = rotation.yaw(or_f[0][0]) @ rotation.pitch(or_f[1][0]) @ rotation.roll(or_f[2][0])
    R_bf = np.linalg.inv(R_fb)
    vel_b = R_bf @ vel_f

    # Rotation from body to aerodynamic frame, Velocity in the aerodynamic frame
    R_ba = rotation.body_aero(vel_b)
    R_ab = np.linalg.inv(R_ba)
    vel_a = R_ab @ vel_b

    #* Force and Torque Calculations
    # Moment Arm from center of mass to center of pressure
    c_p =  2.19 #m
    moment_arm_b = np.array([[c_m - c_p],
                             [0],
                             [0]])
    moment_arm_a = R_ab @ moment_arm_b

    # Density varies with altitude
    rho = atmosphere.density(pos_f[0][0])

    # Total drag coefficient of airframe function imported
    Cd_total = drag_from_csv(pos_f[0][0],vel_b)

    # calculating the sum of aerodynamic forces on the rocket body
    v_mag = np.linalg.norm(vel_a)

    #Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D)

    #Calculate aerodynamic forces on the rocket and the acceleration in the aerodynamic frame
    F_a = -((rho* (v_mag**2) * Sref_a * Cd_total) / 2) * (vel_a/v_mag)
    accel_a = F_a/constants.m0

    # Calculate Torque from Aerodynamic Forces in the Aerodynamic Frame and convert to body frame
    torque_a = np.cross(moment_arm_a, F_a, axis=0)
    torque_b = R_ba @ torque_a

    #* Translational + Angular Acceleration Calculations
    # Calculate angular acceleration of the rocket in the body and fixed frames
    angaccel_b = np.linalg.inv(I) @ torque_b
    angaccel_f = R_fb @ angaccel_b

    # Calculate acceleration in the body frame, convert to fixed frame and calculate total acceleration
    accel_b = R_ba @ accel_a
    accel_f = R_fb @ accel_b - np.array([[constants.g],[0],[0]])

    #* End simulation if rocket reached apogee
    if (np.sign(vel_f[0][0]) == -1):
        break

    #* Updating Values
    # # Append Values to the Arrays

    dic["pos_vals"]["x"].append(float(pos_f[0]))
    dic["pos_vals"]["y"].append(float(pos_f[1]))
    dic["pos_vals"]["z"].append(float(pos_f[2]))

    dic["or_vals"]["Yaw"].append(float(or_f[0]))
    dic["or_vals"]["Pitch"].append(float(or_f[1]))
    dic["or_vals"]["Roll"].append(float(or_f[2]))

    dic["vel_vals"]["Vx"].append(float(vel_f[0]))
    dic["vel_vals"]["Vy"].append(float(vel_f[1]))
    dic["vel_vals"]["Vz"].append(float(vel_f[2]))

    dic["angvel_vals"]["Yaw rate"].append(float(angvel_f[0]))
    dic["angvel_vals"]["Pitch rate"].append(float(angvel_f[1]))
    dic["angvel_vals"]["Roll rate"].append(float(angvel_f[2]))

    dic["accel_vals"]["Ax"].append(float(accel_f[0]))
    dic["accel_vals"]["Ay"].append(float(accel_f[1]))
    dic["accel_vals"]["Az"].append(float(accel_f[2]))

    dic["angaccel_vals"]["Yaw accel"].append(float(angaccel_f[0]))
    dic["angaccel_vals"]["Pitch accel"].append(float(angaccel_f[1]))
    dic["angaccel_vals"]["Roll accel"].append(float(angaccel_f[2]))

    #TODO: Add to dict
    Cd_list.append(Cd_total)
    ref_a_vels.append(Sref_a)

    # Calculate new angular rates and orientation using current values
    or_f = or_f + angvel_f*dt + (0.5 * (angaccel_f * (dt**2)))
    angvel_f = angvel_f + angaccel_f*dt

    # Calculate new velocities and positions using current values
    pos_f = pos_f + (vel_f * dt) + (0.5 * (accel_f * (dt**2)))
    vel_f = vel_f + accel_f*dt

#Print Apogee and total time taken
print("APOGEE (ft):", conversion.m_to_ft(max(dic["pos_vals"]["x"])))
print("Total Time Taken (s):", t)

#* --------------------------------- Plotting --------------------------------- #

#Calculate the number of steps simulated and create a new linspace
simulated_steps = int(total_steps * ((t - start_time) / (end_time - start_time)))
time_flight = np.linspace(start_time,t,simulated_steps,endpoint=False)

plt.plot(time_flight,dic["or_vals"]["Yaw"], label="Yaw Rate", linewidth = 3)
plt.plot(time_flight,dic["or_vals"]["Pitch"], label="Pitch Rate", linewidth = 3)
plt.plot(time_flight,dic["or_vals"]["Roll"], label="Roll Rate")
plt.ylabel("Rad",fontsize = 18); plt.xlabel("Time", fontsize = 18)
plt.xticks(fontsize = 14);plt.yticks(fontsize = 14)
plt.legend(fontsize = 20)
plt.show()

plot.plot_accel_time(dic["accel_vals"], time_flight)


# # #? Coefficient of Drag Plot
# plt.figure(dpi = 200)
# time = np.linspace(0,30,len(Cd_list),endpoint=False)
# plt.plot(time, Cd_list)
# plt.xlabel("Time");plt.ylabel("CD")
# plt.show()

# plot.plot_3d_est(pos_vals, dt, True)
# # Plot yaw
# plt.figure(dpi = 200)
# yaw_vals = []
# accel_vals = []
# for x in np.arange(0,len(or_vals)):
#     yaw_vals.append(or_vals[x][1][0])

# plt.plot(time_flight,yaw_vals)
# plt.ylabel("Yaw (radians)"); plt.xlabel("Time (seconds)")/
# plt.show()
# print(yaw_vals)
# Plot acceleration
# plot.plot_accel_time(accel_vals,time_flight)

# plotting reference area of the rocket against time
# plt.figure(dpi = 200)
# time = np.linspace(0,30,len(ref_a_vels))
# plt.plot(time,ref_a_vels)
# plt.xlabel("Time"); plt.ylabel("Reference area (m^2)")
# plt.show()
