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
D = 0.1524 # meters #* From last launch
# reference area of the rocket 
Sref_a = np.pi*(D/2)**2 


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
current_position = np.array([[],
                           [],
                           []]) # X, Y, Z

current_orientation = np.array([[],
                       [],
                       []]) # Yaw, Pitch, Roll

current_velocity = np.array([[],
                       [],
                       []]) # Vx, Vy, Vz

current_angular_rates = np.array([[],
                       [],
                       []]) # Yaw rate, Pitch rate, Roll rate

#*** Lists for storing values during simulation
h_vals = []
v_vals = []
a_vals = []
omega_vals = []
alpha_vals = []
l1_vals = []
l2_vals = []
power1_vals = []
power2_vals = []

#*** Time setup
time = np.linspace(0,30,10000,endpoint=False)
dt = time[1] - time[0]

#*** Define Max and Min values for Flap Actuation
l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

#* Define initial flap length at start of coasting
l1 = 0
l2 = 0

for t in time:
    #? 1) convert fixed frame velocity to body frame and then to aerodynamic frame using rotational matrix

    #? 2) calculate the aerodynamic forces in the aerodynamic frame using the velocities calculated in 1) 
    #?      Take into account the rotation of the rocket (MAKE A FUNCTION FOR IT!!!): Reference area changes when the rocket rotates
    #?      !!!roll, pitch, yaw are defined in the fixed frame!!!

    #?          a) we have T = I * alpha where alpha is the angular rates (change of roll, pitch, yaw w.r.t time)
    #?             Torque T is calculated from the aerodynamic forces, moment arm is from the center of mass to center of pressure (body frame)
    #?             To do this, convert the moment arm from the body frame to the aerodynamic frame first
    #?             Torque results in a 3x1 matrix in the aerodynamic frame 

    #?          b) Convert this torque to the body frame, now it is a 3x1 matrix where moment of inertia I is a 3x3 matrix (body frame)

    #?          c) invert I and multiply by the new T in the fixed frame to find a 3x1 matrix for the angular rate alpha, this gives alpha in (body frame)
    
    #?          d) convert alpha (body frame) to fixed frame

    #? 3) calculate the acceleration in the aerodynamic frame, this should be just one value 1x1 matrix 

    #? 4) convert 3) to fixed frame to obtain ax, ay, az and propagate x, y, z (m) values or vx, vy, vz values in fixed frame 

    #? 5) propagate roll, pitch, yaw angles using 2)b) 

    #* Convert values to the body frame
    #* Rotation matrix to go from fixed to body frame
    R_fb = rotation.yaw(current_orientation[0][0]) @ rotation.pitch(current_orientation[1][0]) @ rotation.roll(current_orientation[2][0])  
    #* Convert values to the aerodynamic frame and apply forces
    
    #* Recalculate values
    
    

    #* Recalculate Acceleration at each time-step
    a_t1 = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - ((rho*(v_t**2)*Cd_f*W_f*(l1 + l2)) / m) - g

    #* Initial Velocity + Change in Velocity due to instantaneous acceleration
    v_t1 = v_t + (a_t * dt) 

    #* Initial Height + Change in Height from Velocity, Acceleration
    h_t1 = h_t + (v_t * dt) + (0.5 * (a_t *  (dt** 2))) 
    
    #* Roll Rate Calculation
    alpha_t1 = rho * (v_t**2) * Cl_f * W_f * (D*(l1 - l2) + l1**2 - (l2**2))
    
    omega_t1 = omega_t + (alpha_t * dt)

    #* Add values to list
    h_vals.append(h_t)
    v_vals.append(v_t)
    a_vals.append(a_t)
    omega_vals.append(omega_t)
    alpha_vals.append(alpha_t)
    l1_vals.append(l1)
    l2_vals.append(l2)
   
    #* Update new values
    h_t = h_t1
    v_t = v_t1
    a_t = a_t1
    omega_t = omega_t1
    alpha_t = alpha_t1

    #* End simulation if rocket has hit the ground
    if (h_t <= 0):
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