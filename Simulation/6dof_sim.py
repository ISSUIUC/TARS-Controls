import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix

#? import Helper Function Library 
from helper_library import *

# #* EOMS
# hdot = u
# hddot = -((rho * u**2 * Sref_a * Cd_airframe) / (2*m)) - ((rho * u**2 * Cd_flap * W_flap * (l1+l2))/m) - g
# omegadot = (rho * u**2 * Cl_f * W_flap * (l1*D + l1**2 - l2*D - l2**2)) / (2*Ixx)
# f_sym = Matrix.vstack(Matrix([[hdot],[hddot],[omegadot]]))

#* Constants
# Mass of the rocket
m = 21.22  #kg  
# acceleration of gravity
g = 9.81 #m/s^2
#  lift coefficient of the rocket 
Cl_airframe = 0.063 #* From last launch
# lift coefficient of the flaps 
Cl_flap = 2*np.pi*np.sin(45) #*From last launch
# drag coefficient of the rocket 
Cd_airframe = 0.53551
# drag coefficient of the flaps
Cd_flap = 2*np.pi*np.sin(45) #* From last launch
# Width of the flaps 
W_flap = 0.0254 # meters #* From last launch
# Diameter of the rocket 
D = 0.102 # meters #* From last launch
# reference area of the rocket 
Sref_a = np.pi*(D/2)**2 
# Length of the rocket (nosecone + body tube) (m)
l_rocket = 3.02

#TODO: Need to incorporate the density function into the loop while propagating the altitute 
# Density of air 
rho = 1.225 # kg/m^3 #* Changes depending on altitude

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
                  [constants.lateral_velocity],
                  [0]]) # Vx, Vy, Vz

angvel_f = np.array([[np.deg2rad(constants.yaw_rate)],
                     [np.deg2rad(constants.pitch_rate)],
                     [np.deg2rad(constants.roll_rate)]]) # Yaw rate, Pitch rate, Roll rate

# Initialize lists to store values at all time steps
pos_vals = []
or_vals = []
vel_vals =[]
angvel_vals = []
accel_vals = []
angaccel_vals = []

# Time setup
time = np.linspace(0,30,10000,endpoint=False)
# time = np.linspace(0,2,5,endpoint=False)
dt = time[1] - time[0]

# Define max and min values for flap actuation
l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

# Define initial flap length at start of coasting
l1 = 0
l2 = 0

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I,c_m = inertia.I_new(0,0)

for t in time:
    
    #* Velocity Conversions
    # Rotation from fixed to body frame, Velocity converstion to the body frame
    R_fb = rotation.yaw(or_f[0][0]) @ rotation.pitch(or_f[1][0]) @ rotation.roll(or_f[2][0])  
    V_b = R_fb @ vel_f
    
    # Rotation from body to aerodynamic frame, Velocity in the aerodynamic frame
    R_ba = rotation.body_aero(V_b)
    V_a = R_ba @ V_b
    
    #* Force and Torque Calculations
    # Moment Arm from center of mass to center of pressure
    c_p =  2.19 #m
    moment_arm_b = np.array([[c_m - c_p], 
                             [0],
                             [0]]) 
    moment_arm_a = R_ba @ moment_arm_b

    #TODO: Double check this function
    # Calculate the reference area of the rocket
    Sref_a = sref.sref_body(V_b, l_rocket, D)[1]
    
    #TODO: Function for calculating Drag Coefficient based on Reynolds Number - IN PROGRESS
    # Calculate Aerodynamic Forces and Acceleration in Aerodynamic Frame
    F_a = -((rho*np.square(V_a)*Sref_a*Cd_airframe)/2) - (rho*np.square(V_a)*Cd_flap*W_flap*(l1 + l2)) 
    
    accel_a = F_a/m
    
    # Calculate Torque from Aerodynamic Forces in the Aerodynamic Frame and convert to body frame
    torque_a = np.cross(moment_arm_a, F_a, axis=0)
    torque_b = np.linalg.inv(R_ba) @ torque_a    
    

    #* Translational + Angular Acceleration Calculations
    # Calculate angular acceleration of the rocket in the body and fixed frames
    #TODO: Double check this inverse conversion 
    angaccel_b = np.linalg.inv(I) @ torque_b
    angaccel_f = np.linalg.inv(R_fb) @ angaccel_b
 
    # Calculate acceleration in the body frame, convert to fixed frame and calculate total acceleration
    accel_b = np.linalg.inv(R_ba) @ accel_a
    accel_f = np.linalg.inv(R_fb) @ accel_b - np.array([[g],[0],[0]])
    
    
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
    
    # Calculate new angular rates and orientation using current values
    or_f = or_f + angvel_f*dt + (0.5 * (angaccel_f * (dt**2)))
    angvel_f = angvel_f + angaccel_f*dt

    # Calculate new velocities and positions using current values
    pos_f = pos_f + (vel_f * dt) + (0.5 * (accel_f * (dt**2)))
    vel_f = vel_f + accel_f*dt

