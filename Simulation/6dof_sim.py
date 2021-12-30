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
start_alt = conversion.ft_to_m(727.65)
start_vel = conversion.ft_to_m(491.29)
v_t = start_vel
h_t = start_alt
a_t = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - g
omega_t = 0.5
alpha_t = 0.1

start_position = np.array([[],
                           [],
                           []])

start_orientation = np.array([[],
                       [],
                       []])

start_velocity = np.array([[],
                       [],
                       []])

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

#* Plots of states and inputs
fig, axs = plt.subplots(2,2,figsize=(10,5))
axs[0,0].plot(time[0:len(h_vals)],h_vals,label="Altitude")
axs[0,0].set_ylabel("Alt + VV (m, m/s)")
axs[0,0].plot(time[0:len(h_vals)],v_vals,label="Flap Length",color="orange")
axs[0,1].yaxis.set_label_position("right")
axs[0,1].yaxis.tick_right()
axs[0,1].set_ylabel("Flap Length (m)")
axs[0,1].plot(time[0:len(h_vals)], l1_vals, color="tab:purple")
axs[0,1].plot(time[0:len(h_vals)], l2_vals, color="tab:pink")
axs[1,0].plot(time[0:len(h_vals)],omega_vals,label="Roll Rate",color="tab:red")
axs[1,0].set_ylabel("Angular Velocity (rad/s)")
axs[1,1].plot(time[0:len(h_vals)],alpha_vals,label="Roll Acceleration",color="tab:green")
axs[1,1].yaxis.set_label_position("right")
axs[1,1].yaxis.tick_right()
axs[1,1].set_ylabel("Angular Acceleration (rad/s^2)")
axs[1,0].set(xlabel="Time (seconds)")
axs[1,1].set(xlabel="Time (seconds)")
plt.show()