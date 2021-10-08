import numpy as np
import sympy as sym
import json
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.interpolate import interp1d
from sympy import Matrix

def lqr(A, B, Q, R):
    P = linalg.solve_continuous_are(A, B, Q, R)
    K = linalg.inv(R) @  B.T @ P
    return K

def ft_to_m(measurement):
    return (measurement / 3.2808) 

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
Q_diag = np.array([0,0,1]) * 1/100  #h, hdot, omega
Q = np.diag(Q_diag)

## Designing a R matrix
R_diag =  10000 * np.array([1,1]) #l1, l2
R = np.diag(R_diag)

# K matrix
K = lqr(A,B,Q,R)
print(K)

# K = np.array( [[0,0,5],
#               [0,0,-5]]) * 5/10000

# print(K)

## testing to see if the A - BK has all negative eigenvalues 
F = A - B@K
eig = np.linalg.eigvals(F)
print(eig.real < 0)

#* Sim
#*** Initial State at Burnout
start_alt = ft_to_m(727.65)
start_vel = ft_to_m(491.29)
v_t = start_vel
h_t = start_alt
a_t = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - g
omega_t = -0.5
alpha_t = -0.1

#*** Lists for storing values during simulation
h_vals = []
v_vals = []
a_vals = []
omega_vals = []
alpha_vals = []
l1_vals = []
l2_vals = []

#*** Time setup
time = np.linspace(0,30,10000,endpoint=False)
dt = time[1] - time[0]

#*** Define Max and Min values for Flap Actuation
l_max = ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation, 0 inches is minimum

#* Starting Actuation Length
l1 = 0
l2 = 0

for t in time:

    #* Recalculate Acceleration at each time-step
    a_t1 = -((rho*(v_t**2)*Sref_a*Cd_a) / (2*m)) - ((rho*(v_t**2)*Cd_f*W_f*(l1 + l2)) / m) - g

    #* Initial Velocity + Change in Velocity due instantaneous acceleration
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

    #* Control Input Calculation
    x = [[h_t],[v_t],[omega_t]]
    u = -K @ x
    l1 = u[0][0]
    l2 = u[1][0]

    #* Control Input Damping
    if (l1 > l_max):
        l1 = l_max
    elif (l1 < l_min):
        l1 = l_min

    if (l2 > l_max):
        l2 = l_max
    elif (l2 < l_min):
        l2 = l_min

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
axs[0,0].set_ylabel("Alt + Vert Vel (m, m/s)")
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