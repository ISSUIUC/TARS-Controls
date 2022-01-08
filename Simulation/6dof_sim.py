import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix

#? import Helper Function Library 
from helper_library import *

# #* EOMS
# hdot = u
# hddot = -((rho * u**2 * Sref_a * Cd_a) / (2*m)) - ((rho * u**2 * Cd_f * W_f * (l1+l2))/m) - g
# omegadot = (rho * u**2 * Cl_f * W_f * (l1*D + l1**2 - l2*D - l2**2)) / (2*Ixx)
# f_sym = Matrix.vstack(Matrix([[hdot],[hddot],[omegadot]]))

#* Constants
# Mass of the rocket
m = 21.22  #kg  
# moment of inertia
Ixx = 0.5*m*0.1524**2 #* From last launch
# acceleration of gravity
g = 9.81 #m/s^2
#  lift coefficient of the rocket 
Cl_a = 0.063 #* From last launch
# lift coefficient of the flaps 
Cl_f = 2*np.pi*np.sin(45) #*From last launch
# drag coefficient of the rocket 
Cd_a = 0.70649 
# drag coefficient of the flaps
Cd_f = 2*np.pi*np.sin(45) #* From last launch
# Width of the flaps 
W_f = 0.0254 # meters #* From last launch
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
##* Fixed Frame(F) - fixed to the launch rail, not taking into account rotation of the Earth
# X points vertically
# Y is chosen through the right-hand rule
# Z points towards the launch rail


##* Body-Fixed Frame(B) - Origin at the center of mass of the rocket
# X points through the nose of the rocket
# Y is chosen through the right-hand rule
# Z points towards the launch lugs on the rocket


##* Aerodynamic Frame (Body-Fixed) with Origin at the center-of-mass
# X aligned with aerodynamic velocity
# Y is chosen through the right-hand rule
# Z lies in plane that goes through the centerline of the rocket and the launch lug

#* Sim
#*** Initial State at Burnout
# start_alt = conversion.ft_to_m(727.65)
# start_vel = conversion.ft_to_m(491.29)
# v_t = start_vel
# h_t = start_alt
# a_t = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - g
# omega_t = 0.5
# alpha_t = 0.1

#* NumPy arrays to store current states
#? Fixed
current_position = np.array([[],
                           [],
                           []]) # X, Y, Z
#? Fixed
current_orientation = np.array([[],
                       [],
                       []]) # Yaw, Pitch, Roll
#? Fixed
current_velocity = np.array([[],
                       [],
                       []]) # Vx, Vy, Vz
#? Fixed
current_angular_rates = np.array([[],
                       [],
                       []]) # Yaw rate, Pitch rate, Roll rate
#? Fixed
current_acceleration =  np.array([[],
                       [],
                       []]) # ax, ay, az
#? Fixed
current_alpha =  np.array([[],
                       [],
                       []]) # change of angular rate w.r.t time

#*** Lists for storing values during simulation
h_vals = []
v_vals = []
a_vals = []
orientation_vals = []
ang_vel_vals = []
alpha_vals = []



#*** Time setup
time = np.linspace(0,30,10000,endpoint=False)
dt = time[1] - time[0]

#*** Define Max and Min values for Flap Actuation
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
    R_fb = rotation.yaw(current_orientation[0][0]) @ rotation.pitch(current_orientation[1][0]) @ rotation.roll(current_orientation[2][0])  
    V_b = R_fb @ current_velocity
    
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
    moment_arm = np.array([[c_m - c_p], 
                           [0],
                           [0]]) # Body Frame
    moment_arm_aero = R_ba @ moment_arm
        
    #TODO: Double check this function

    Sref_a = sref.sref_body(V_b, D, l_rocket)
    
    #TODO: Function for calculating Drag Coefficient based on Reynolds Number - IN PROGRESS
    
    #* Calculate Forces in Aerodynamic Frame
    F_aero = -((rho*np.square(V_a)*Sref_a*Cd_a)/2) - (rho*np.square(V_a)*Cd_f*W_f*(l1 + l2))
    #* Calculate Acceleration in Aerodynamic Frame (1X1)
    a_aero = F_aero/m
    #* Torque from Aerodynamic Forces in the Aerodynamic Frame
    torque_aero = np.cross(moment_arm_aero, F_aero)
    
    #* DONT FORGET GRAVITY

    #?    b) Convert this torque to the body frame, now it is a 3x1 matrix where moment of inertia I is a 3x3 matrix (body frame)

    torque_body = np.linalg.inv(R_ba) @ torque_aero

    #?    c) invert I and multiply by the new T in the fixed frame to find a 3x1 matrix for the angular acceleration alpha, this gives alpha in (body frame)
    
    alpha_body = np.linalg.inv(I) @ torque_body
    
    #?    d) convert alpha (body frame) to fixed frame

    #TODO: Double check this inverse conversion 
    alpha_fixed = np.linalg.inv(R_fb) @ alpha_body

    #? 3) convert a_aero to fixed frame to obtain ax, ay, az and propagate x, y, z (m) values or vx, vy, vz values in fixed frame 

    a_body = np.linalg.inv(R_ba) @ a_aero
    a_fixed = np.linalg.inv(R_fb) @ a_body - np.array([[0],[0],[g]])
    
    #? 4) propagate roll, pitch, yaw angles using alpha_fixed

    new_angular_rates = current_angular_rates + alpha_fixed*dt
    new_orientation = current_orientation + current_angular_rates*dt + (0.5 * (alpha_fixed * (dt**2)))

    #? 5) propagate vx, vy, vz

    new_velocity = current_acceleration + current_velocity*dt

    #? 5) propagate x, y, z

    new_position = current_position + (current_velocity * dt) + (0.5 * (current_acceleration * (dt**2)))

    #? 6) Add values to list 

    h_vals.append(new_position)
    v_vals.append(new_velocity)
    a_vals.append(a_fixed)
    ang_vel_vals.append(new_angular_rates)
    orientation_vals.append(new_orientation)
    alpha_vals.append(alpha_fixed)

   
    #* Update new values
    current_velocity = new_velocity
    current_position = new_position
    current_acceleration = a_fixed
    current_angular_rates = new_angular_rates
    current_alpha = alpha_fixed
    current_orientation = new_orientation

    #* End simulation if rocket reached apogee
    if (np.sign(new_velocity[0][2]) == -1):
        break

# #* Plots of states and inputs
# fig, axs = plt.subplots(2,2,figsize=(10,5))
# axs[0,0].plot(time[0:len(h_vals)],h_vals,label="Altitude")
# axs[0,0].set_ylabel("Alt + VV (m, m/s)")
# axs[0,0].plot(time[0:len(h_vals)],v_vals,label="Flap Length",color="orange")
# axs[0,1].yaxis.set_label_position("right")
# axs[0,1].yaxis.tick_right()
# axs[0,1].set_ylabel("Flap Length (m)")
# axs[0,1].plot(time[0:len(h_vals)], l1_vals, color="tab:purple")
# axs[0,1].plot(time[0:len(h_vals)], l2_vals, color="tab:pink")
# axs[1,0].plot(time[0:len(h_vals)],omega_vals,label="Roll Rate",color="tab:red")
# axs[1,0].set_ylabel("Angular Velocity (rad/s)")
# axs[1,1].plot(time[0:len(h_vals)],alpha_vals,label="Roll Acceleration",color="tab:green")
# axs[1,1].yaxis.set_label_position("right")
# axs[1,1].yaxis.tick_right()
# axs[1,1].set_ylabel("Angular Acceleration (rad/s^2)")
# axs[1,0].set(xlabel="Time (seconds)")
# axs[1,1].set(xlabel="Time (seconds)")
# plt.show()