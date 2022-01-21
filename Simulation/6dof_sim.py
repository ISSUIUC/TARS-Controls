from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
import pandas as pd

#* Import Helper Function Library 
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation

#* Importing RasAero Packege
rasaero = pd.read_csv("Simulation/Src/RASAero.csv")
# extracting the columns of interest 
mach_num = rasaero.mach.values; aoa = rasaero.alpha_deg.values; cd = rasaero.cd_power_off.values; protub = rasaero.protuberance.values
# narrowing down the columns using mach number range (0.02 - 1.01)
min_index = min(np.where(mach_num == 0.01)[0])
max_index = min(np.where(mach_num == 1.01)[0])
# re-make the lists using this range 
mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]

# Python function to print common elements in three sorted arrays

def drag_from_csv(z, velocity_body, pro):
    mach = round(np.linalg.norm(velocity_body) / atmosphere.speed_sound(z),2)
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[1][0]
    vz_b = velocity_body[2][0]
    alpha = abs(round(np.rad2deg(np.arctan2(vz_b,vx_b)), 0))

    print(mach,alpha)


    mach_index_array = np.where(mach_num == mach)[0]; alpha_index_array = np.where(aoa == alpha)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))

    if len(intersection) == 0:
        Cd_list.append(Cd_list[-1])
        return Cd_list[-1]

    index = intersection[0]
    Cd_list.append(cd[index])
    return cd[index]
 

#TODO: Debugging
#Check initial acceleration from simulation matches up with OpenRocket/RASAero #! Doesnt look like it
#Do we need different reference areas for drag? Should it be a vector?

#* Constants
# Mass of the rocket (dry)
m = 21.22  #kg  
# acceleration of gravity
g = 9.81 #m/s^2
# length of rocket 
l_rocket = 3.02

# nosecone angle (rad)
angle = 0.069189
###! Ignore this if this doesn't work 
#* Total length of Rocket 
l = 3.0226
#* Rocket outer diameter 
D = 0.1056132
d_b = D
d_d = D
#* Body Tube Length 
L_b = 2.2352
#* Nose Cone Length 
L_n = 0.762
#* Fin thickness
T_f = 0.0029972
#* true length of the fin from inner to outer edge/ root chord 
L_m = 0.2032
#* Number of fins 
n = 3
#* fin platform area 
A_fp = 0.011532235
#* fin height 
d_f = 0.08255
###! Ignore if this doesn't work

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
pos_f = np.array([[constants.x],
                  [0],
                  [0]]) # X, Y, Z

or_f = np.array([[0],
                 [0],
                 [0]]) # Yaw, Pitch, Roll

vel_f = np.array([[constants.vx],
                  [0],
                  [0]]) # Vx, Vy, Vz

angvel_f = np.array([[constants.yaw_rate],
                     [constants.pitch_rate],
                     [constants.roll_rate]]) # Yaw rate, Pitch rate, Roll rate

# Initialize lists to store values at all time steps
pos_vals = []
or_vals = []
vel_vals =[]
angvel_vals = []
accel_vals = []
angaccel_vals = []
ref_a_vels = []

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
I,c_m = rocket.I_new(0,0)
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

    #TODO: Double check this function
    # Calculate the reference area of the rocket
    Sref_a = rocket.sref_body(vel_b, l_rocket, D)[1]
    
    #TODO: Function for calculating Drag Coefficient based on Reynolds Number - IN PROGRESS
    
    # Varying density function imported 
    rho = atmosphere.density(pos_f[0][0])
    
    # Total drag coefficient of airframe function imported 
    # Cd_total = rocket.total_drag_scaled(pos_f[0][0],l_rocket,D,vel_b, Sref_a, angle)
    Cd_total = drag_from_csv(pos_f[0][0],vel_b,0)
    # inter = drag_from_csv(pos_f[0][0],vel_b,0)
    
    # print(inter)
    # print(mach,alpha, pro)
    Cd_list.append(Cd_total)
    
    # calculating the sum of aerodynamic forces on the rocket body
    v_mag = np.linalg.norm(vel_a)
    F_a = -((rho*(v_mag**2)*Sref_a*Cd_total)/2)*(vel_a/(np.linalg.norm(vel_a)))
    accel_a = F_a/m

    # Calculate Torque from Aerodynamic Forces in the Aerodynamic Frame and convert to body frame
    torque_a = np.cross(moment_arm_a, F_a, axis=0)
    torque_b = R_ba @ torque_a    
    
    #* Translational + Angular Acceleration Calculations
    # Calculate angular acceleration of the rocket in the body and fixed frames
    #TODO: Double check this inverse conversion 
    angaccel_b = np.linalg.inv(I) @ torque_b
    angaccel_f = R_fb @ angaccel_b
 
    # Calculate acceleration in the body frame, convert to fixed frame and calculate total acceleration
    accel_b = R_ba @ accel_a
    accel_f = R_fb @ accel_b - np.array([[g],[0],[0]])
    
    #* End simulation if rocket reached apogee
    if (np.sign(vel_f[0][0]) == -1):
        break
    
    #* Updating Values
    # # Append Values to the Arrays
    pos_vals.append(pos_f)
    or_vals.append(or_f)
    vel_vals.append(vel_f)
    angvel_vals.append(angvel_f)
    accel_vals.append(accel_f)
    angaccel_vals.append(angaccel_f)
    Cd_list.append(Cd_total)
    ref_a_vels.append(Sref_a)


    
    # Calculate new angular rates and orientation using current values
    or_f = or_f + angvel_f*dt + (0.5 * (angaccel_f * (dt**2)))
    angvel_f = angvel_f + angaccel_f*dt

    # Calculate new velocities and positions using current values
    pos_f = pos_f + (vel_f * dt) + (0.5 * (accel_f * (dt**2)))
    vel_f = vel_f + accel_f*dt

print(max(pos_vals[-1][0]))

# plotting reference area of the rocket against time
# plt.figure(dpi = 200)
# time = np.linspace(0,30,len(ref_a_vels))
# plt.plot(time,ref_a_vels)
# plt.xlabel("Time"); plt.ylabel("Reference area (m^2)")
# plt.show()


# # # checking the yaw pitch roll values
yaw_vals = []
pitch_vals = []
roll_vals = []
pitch_rate = []
yaw_rate = []
roll_rate = []
time = np.linspace(0,30,len(or_vals),endpoint=False)
for x in np.arange(0,len(angvel_vals)):
    pitch_vals.append(or_vals[x][1][0])
    yaw_vals.append(or_vals[x][0][0])
    roll_vals.append(or_vals[x][2][0])
    yaw_rate.append(angvel_vals[x][0][0])
    pitch_rate.append(angvel_vals[x][1][0])
    roll_rate.append(angvel_vals[x][2][0])


plt.plot(time,yaw_vals, label="Yaw Rate"); plt.plot(time,pitch_vals, label="Pitch Rate"); #plt.plot(time,roll_vals, label="Roll Rate")
plt.ylabel("Rad/s"); plt.xlabel("Time")
plt.legend()
plt.show()

# #? Coefficient of Drag Plot 
plt.figure(dpi = 200)
time = np.linspace(0,30,len(Cd_list),endpoint=False)
plt.plot(time, Cd_list)
plt.xlabel("Time");plt.ylabel("CD")
plt.show()


#Calculate the number of steps simulated before break 
# simulated_steps = int(total_steps * ((t - start_time) / (end_time - start_time))) 
# time_flight = np.linspace(start_time,t,simulated_steps,endpoint=False)

# plot.plot_3d_est(pos_vals, dt, True)
# # # Plot yaw
# plt.figure(dpi = 200) 
# yaw_vals = []
# accel_vals = []
# for x in np.arange(0,len(or_vals)):
#     yaw_vals.append(or_vals[x][1][0])
    
# plt.plot(time_flight,yaw_vals)
# plt.ylabel("Yaw (radians)"); plt.xlabel("Time (seconds)")
# plt.show()
# print(yaw_vals)
# Plot acceleration
# plot.plot_accel_time(accel_vals,time_flight)


