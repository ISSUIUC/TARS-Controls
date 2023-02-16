import numpy as np
import sympy as sym
import json
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.interpolate import interp1d
from sympy import Matrix
import control
from statistics import mean
import random
from scipy import stats

def lqr(A, B, Q, R):
    P = linalg.solve_continuous_are(A, B, Q, R)
    K = linalg.inv(R) @  B.T @ P
    return K

def ft_to_m(measurement):
    return (measurement / 3.2808) 


def k_matrix(x1,x2,x3):
    return np.array([[x1,x2,x3],[-x1,-x2,-x3]])

# defining components of position 
h, y, z = sym.symbols('h, y, z')
# defining roll, pitch and yaw
psi, theta, phi  = sym.symbols('psi, theta, phi')
# defining the components of the linear velocity 
u, v, w = sym.symbols('u, v, w')
# defining compoenents for angular velocity 
wx, wy, wz = sym.symbols('wx, wy, wz')
# defining inputs
l1, l2 = sym.symbols('l1, l2')
# defining the physical parameters on the rocket 
g, W_f, Ixx, Sref_a, rho, m, D = sym.symbols('g, W_f, Ixx, Sref_a, rho, m, D')
# defining the coefficients 
Cd_a, Cd_f, Cl_a, Cl_f = sym.symbols('Cd_a, Cd_f, Cl_a, Cl_f')

#* EOMS
hdot = u
hddot = -((rho * u**2 * Sref_a * Cd_a) / (2*m)) - ((rho * u**2 * Cd_f * W_f * (l1+l2))/m) - g
omegadot = (rho * u**2 * Cl_f * W_f * (l1*D + l1**2 - l2*D - l2**2)) / (2*Ixx)
f_sym = Matrix.vstack(Matrix([[hdot],[hddot],[omegadot]]))

#* State Space Model
s = [h, hdot, wx]
i = [l1, l2]
p = [Cd_f, Sref_a,Cd_a,rho, W_f, m, g, Cl_f, Ixx, D]
f = sym.lambdify(s + i + p, f_sym)

#* Plug in Constants
# Mass of the rocket
m = 50 #kg
# moment of inertia
Ixx = 0.5*m*0.1524**2
# acceleration of gravity
g = 9.81 #m/s^2
#  lift coefficient of the rocket 
Cl_a = 0.063
# lift coefficient of the flaps 
Cl_f = 2*np.pi*np.sin(45)
# drag coefficient of the rocket 
Cd_a = 0.58
# drag coefficient of the flaps
Cd_f = 2*np.pi*np.sin(45)
# Width of the flaps 
W_f = 0.0254 # meters
# Diameter of the rocket 
D = 0.1524 # meters
# reference area of the rocket 
Sref_a = np.pi*(D/2)**2
# Density of air 
rho = 1.225 # kg/m^3

#* Define and Plug in Equilibrium Values
h_e = ft_to_m(900) 
hdot_e = ft_to_m(400) 
omega_e = 0
l1_e = ft_to_m(1/12)
l2_e = ft_to_m(1/12)

s_eq = [h_e,hdot_e,omega_e]
i_eq = [l1_e,l2_e]
p_eq = [Cd_f,Sref_a,Cd_a,rho,W_f,m,g,Cl_f,Ixx,D]

feq = f(*s_eq,*i_eq,*p_eq)
# print("XDot At Equilbrium: ")
# print(feq)

#* Compute A and B
A_sym = f_sym.jacobian(s)
A_num = sym.lambdify(s + i + p,A_sym)
A = A_num(*s_eq, *i_eq, *p_eq)

B_sym = f_sym.jacobian(i)
B_num = sym.lambdify(s + i + p,B_sym)
B = B_num(*s_eq, *i_eq, *p_eq)

#* Compute Controllability
n = len(B)

# Initialize the Controllability Matrix (W)
W = B

# Create the Controllability Matrix (W)
for i in range(1,n):
    new_mat = np.linalg.matrix_power(A, i) @ B
    W = np.block([W,new_mat])

# Make sure that the rank of the matrix is equal to the number of states
# if (np.linalg.matrix_rank(W) == n):
#     print("Full rank")
# else:
#     print("Rank deficient") 

#* LQR Design
## Designing a Q matrix 
Q_diag = np.array([0,0,100]) #h, hdot, omega
Q = np.diag(Q_diag)

## Designing a R matrix
R_diag =  3 * np.array([1,1]) #l1, l2
R = np.diag(R_diag)

# K matrix
K = lqr(A,B,Q,R)

# print(K)

# K = [[0,0,1],
#      [0,0,1]]

K = np.array([[0,0,0.005],
              [0,0,-0.005]])

#? Defining a set of negative eigenvalues 
# p = [0,-1000,0]
# K = control.place(A,B,p);
# print(K)

## testing to see if the A - BK has all negative eigenvalues
print(B)

F = A - B@K
eig = np.linalg.eigvals(F)
print(eig)
print(eig.real < 0)

#* Sim
#*** Initial State at Burnout
start_alt = ft_to_m(727.65)
start_vel = ft_to_m(491.29)
v_t = start_vel
h_t = start_alt
a_t = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - g
omega_t = 0.5
alpha_t = 0.1

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
l_max = ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

#* Define initial flap length at start of coasting
l1 = 0
l2 = 0

#* Define a constant torque output from the servo motor 
T = 2.353596 # newton meter 

#* Define power output from the servo (operating at 0.009 amp and 6 V)
power1 = 0.009*6 # watt
power2 = 0.009*6 # watt

#* Generating a random angular velocity as a disturbance (wind)
omega_rand = random.uniform(-0.5,0.5)

#* Goal
omega_goal = 0
h_goal = 3000
v_goal = 0


k_p = k_matrix(0, 0, 0.003)
k_I = k_matrix(0, 0, 0.0000007)
k_d = k_matrix(0, 0, 0.0000007)
for t in time:
    e_sum = 0

    #* Recalculate Acceleration at each time-step
    a_t1 = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - ((rho*(v_t**2)*Cd_f*W_f*(l1 + l2)) / m) - g

    #* Initial Velocity + Change in Velocity due to instantaneous acceleration
    v_t1 = v_t + (a_t * dt) 

    #* Initial Height + Change in Height from Velocity, Acceleration
    h_t1 = h_t + (v_t * dt) + (0.5 * (a_t *  (dt** 2))) 

     
    #* Roll Rate Calculation
    alpha_t1 = rho * (v_t**2) * Cl_f * W_f * (D*(l1 - l2) + l1**2 - (l2**2))
    
    omega_t1 = omega_t + (alpha_t * dt)

    #* Putting this wind disturbance at a randomly generated time t within the specified range "time"
    for x in np.arange(0,5):
        if t == random.choice(time):
            omega_t1 = omega_rand + (alpha_t * dt)

    #* Add values to list
    h_vals.append(h_t)
    v_vals.append(v_t)
    a_vals.append(a_t)
    omega_vals.append(omega_t)
    alpha_vals.append(alpha_t)
    l1_vals.append(l1)
    l2_vals.append(l2)

    #* Control Input Calculation
    #* Defining error from the timestep = t
    e_prev = np.array([[0],[0],[omega_goal - omega_t]])
    #* Defining error from timestep = t + 1
    e =  np.array([[0],[0],[omega_goal - omega_t1]])
    #* Summing up all the errors 
    e_sum = e_sum + e*dt
    #* Derivative Term (change of the error over time)
    dedt = (e - e_prev)/dt 

    #* PID controller input 
    u = k_p @ e + k_I @ e_sum + k_d @ dedt
    l1_t1 = u[0][0]
    l2_t1 = u[1][0]
  


    e_prev = e

    # x = [[h_t],[v_t],[omega_t]]
    # u = -K @ x
    # l1_t1 = u[0][0]
    # l2_t1 = u[1][0]

    #* Control Input Damping
    if (l1_t1 > l_max):
        l1_t1 = l_max
    elif (l1_t1 < l_min):
        l1_t1 = l_min

    if (l2_t1 > l_max):
        l2_t1 = l_max
    elif (l2_t1 < l_min):
        l2_t1 = l_min
    
    
    #* Calculating the flap extension velocity 
    v_flap1 = abs((l1_t1 - l1)/dt)
    v_flap2 = abs((l2_t1 - l2)/dt)
   
    #* Translating to gear velocites m/s (this is very randomo right now)
    v_g1 = v_flap1*2
    v_g2 = v_flap2*1.5

    #* We need to check whether this pre-defined torque multiplied by 
    #* the velocity will exceed tha maximum power limit
    gear_r = 0.05 # gear redius in meters
    gear_c = 2*np.pi*gear_r # gear circumference

    #* obtaining revolution per second
    gear1_rps = v_g1/gear_c; gear2_rps = v_g2/gear_c 
    
    #* converting to rpm
    gear1_rpm = gear1_rps*60; gear2_rpm = gear2_rps*60

    #* Power Calculation 
    power1_t1 = T*np.pi*gear1_rpm/30
    power2_t1 = T*np.pi*gear2_rpm/30

    power1_vals.append(power1_t1)
    power2_vals.append(power2_t1)
    # print(power2_t1)
    # if power1_t1 > power1 or power2_t1 > power2:
    #     print("New Controller")
    #     break
    #* Update new values
    h_t = h_t1
    v_t = v_t1
    a_t = a_t1
    omega_t = omega_t1
    alpha_t = alpha_t1
    l1 = l1_t1
    l2 = l2_t1


    

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
axs[0,1].plot(time[0:len(h_vals)], l1_vals, color="tab:purple", label="L1 Length");
axs[0,1].plot(time[0:len(h_vals)], l2_vals, color="tab:pink",label="L2 Length"); 
axs[1,0].plot(time[0:len(h_vals)],omega_vals,label="Roll Rate",color="tab:red"); 
axs[1,0].set_ylabel("Angular Velocity (rad/s)")
axs[1,1].plot(time[0:len(h_vals)],alpha_vals,label="Roll Acceleration",color="tab:green")
axs[1,1].yaxis.set_label_position("right")
axs[1,1].yaxis.tick_right()
axs[1,1].set_ylabel("Angular Acceleration (rad/s^2)")
axs[1,0].set(xlabel="Time (seconds)")
axs[1,1].set(xlabel="Time (seconds)")

plt.show()

power1_mean = mean(power1_vals)
power2_mean = mean(power2_vals)
print(f"Average Power Required for servo 1 = {power1_mean} watts")
print(f"Average Power Required for servo 2 = {power2_mean} watts")
print(max(h_vals))
if power2_mean > power2 or power1_mean > power1:
    print("Adjust the controller")