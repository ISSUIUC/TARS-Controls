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
# moment of inertia
# Ixx = 0.5*m*0.1524**2 #* From last launch
# acceleration of gravity
g = 9.81 #m/s^2
#  lift coefficient of the rocket 
Cl_airframe = 0.063 #* From last launch
# lift coefficient of the flaps 
Cl_flap = 2*np.pi*np.sin(45) #*From last launch
# drag coefficient of the rocket 
Cd_airframe = 0.70649 
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

# ---------------------------- Frames we are using --------------------------- #
##* Fixed Frame - fixed to the launch rail, not taking into account rotation of the Earth
# X points vertically
# Y is chosen through the right-hand rule
# Z points towards the launch rail


##* Body-Fixed Frame - Origin at the center of mass of the rocket
# X points through the nose of the rocket
# Y is chosen through the right-hand rule
# Z points towards the launch lugs on the rocket


##* Aerodynamic Frame (Body-Fixed) with Origin at the center-of-mass
# X aligned with aerodynamic velocity
# Y is chosen through the right-hand rule
# Z lies in plane that goes through the centerline of the rocket and the launch lug

# ---------------------------- Naming Conventions ---------------------------- #
#* Values
# pos = position(x,y,z)
# or = orientation(yaw,pitch,roll)
# vel = velocity(x,y,z)
# angvel = angular velocity(yaw,pitch,roll)
# accel = acceleration(x,y,z)
# angaccel = angular acceleration (yaw,pitch,roll)

#* Frames
#f = fixed frame
#b = body frame
#a = aerodynamic frame

#* Time Step
#0 = current time step
#1 = next time step

# pos_f_0 -> current position in the fixed frame
# angvel_b_1 -> next angular velocity in the body frame
#* ------------------------------ Simulation Code ----------------------------- #

#* NumPy arrays to store current states

pos_f_0 = np.array([[],
                    [],
                    []]) # X, Y, Z

or_f_0 = np.array([[],
                   [],
                   []]) # Yaw, Pitch, Roll

vel_f_0 = np.array([[],
                    [],
                    []]) # Vx, Vy, Vz

angvel_f_0 = np.array([[],
                       [],
                       []]) # Yaw rate, Pitch rate, Roll rate

accel_f_0 =  np.array([[],
                       [],
                       []]) # ax, ay, az

angaccel_f_0 =  np.array([[],
                          [],
                          []]) # change of angular rate w.r.t time

#* Lists to store values at all time steps
pos_vals = []
or_vals = []
vel_vals =[]
angvel_vals = []
accel_vals = []
angaccel_vals = []


#* Time setup
time = np.linspace(0,30,10000,endpoint=False)
dt = time[1] - time[0]

#* Define Max and Min values for Flap Actuation
l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

#* Define initial flap length at start of coasting
l1 = 0
l2 = 0

#* Calculate moments of Inertia, center of mass
I,c_m = inertia.I_new(0,0)

for t in time:
    #? 1) convert fixed frame velocity to body frame and then to aerodynamic frame using rotational matrix
    #* Rotation from fixed to body frame, Velocity calculation in the body frame
    R_fb = rotation.yaw(or_f_0[0][0]) @ rotation.pitch(or_f_0[1][0]) @ rotation.roll(or_f_0[2][0])  
    V_b = R_fb @ vel_f_0
    
    #* Rotation from body to aerodynamic frame
    R_ba = rotation.body_aero(V_b)
    V_a = R_ba @ V_b
    
    #? 2) calculate the aerodynamic forces in the aerodynamic frame using the velocities calculated in 1) 
    #?      Take into account the rotation of the rocket (MAKE A FUNCTION FOR IT!!!): Reference area changes when the rocket rotates
    #?      !!!roll, pitch, yaw are defined in the fixed frame!!!

    #?   a) we have T = I * alpha where alpha is the angular acceleration (change of rate of roll, pitch, yaw w.r.t time)
    #?   Torque T is calculated from the aerodynamic forces, moment arm is from the center of mass to center of pressure (body frame)
    #?   To do this, convert the moment arm from the body frame to the aerodynamic frame first
    #?   Torque results in a 3x1 matrix in the aerodynamic frame 
    
    #* Moment Arm from center of mass to center of pressure
    c_p =  2.19 #m
    moment_arm_b = np.array([[c_m - c_p], 
                           [0],
                           [0]]) 
    moment_arm_a = R_ba @ moment_arm_b
        
    #TODO: Double check this function

    Sref_a = sref.sref_body(V_b, D, l_rocket)
    
    #TODO: Function for calculating Drag Coefficient based on Reynolds Number - IN PROGRESS
    
    #* Calculate Aerodynamic Forces in Aerodynamic Frame
    F_a = -((rho*np.square(V_a)*Sref_a*Cd_airframe)/2) - (rho*np.square(V_a)*Cd_flap*W_flap*(l1 + l2))
    #* Calculate Acceleration from Aerodynamic Forces in Aerodynamic Frame 
    accel_a = F_a/m
    
    #* Calculate Torque from Aerodynamic Forces in the Aerodynamic Frame and convert to body frame
    torque_a = np.cross(moment_arm_a, F_a)
    torque_b = np.linalg.inv(R_ba) @ torque_a    

    #* Convert torque from aerodynamic forces into the body frame and calculate angular acceleration of the rocket in the body frame
    
    angaccel_b = np.linalg.inv(I) @ torque_b
    
    #?    d) convert alpha (body frame) to fixed frame

    #TODO: Double check this inverse conversion 
    angaccel_f = np.linalg.inv(R_fb) @ angaccel_b

    #? 3) convert a_aero to fixed frame to obtain ax, ay, az and propagate x, y, z (m) values or vx, vy, vz values in fixed frame 

    accel_b = np.linalg.inv(R_ba) @ accel_a
    accel_f = np.linalg.inv(R_fb) @ accel_body - np.array([[0],[0],[g]])
    
    #? 4) propagate roll, pitch, yaw angles using alpha_f

    new_angular_rates = current_angular_rates + alpha_f*dt
    new_orientation = or_f_0 + current_angular_rates*dt + (0.5 * (alpha_f * (dt**2)))

    #? 5) propagate vx, vy, vz

    new_velocity = current_acceleration + vel_f_0*dt

    #? 5) propagate x, y, z

    new_position = current_position + (vel_f_0 * dt) + (0.5 * (current_acceleration * (dt**2)))

    #? 6) Add values to list 

    # h_vals.append(new_position)
    # v_vals.append(new_velocity)
    # a_vals.append(a_fixed)
    # ang_vel_vals.append(new_angular_rates)
    # orientation_vals.append(new_orientation)
    # alpha_vals.append(alpha_f)

   
    #* Update new values
    vel_f_0 = new_velocity
    current_position = new_position
    current_acceleration = a_fixed
    current_angular_rates = new_angular_rates
    current_alpha = alpha_f
    or_f_0 = new_orientation

    #* End simulation if rocket reached apogee
    if (np.sign(new_velocity[0][2]) == -1):
        break