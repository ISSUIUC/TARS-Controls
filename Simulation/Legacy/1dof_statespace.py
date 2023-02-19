import time
import matplotlib.pyplot as plt
from scipy import linalg
import sympy as sym
import numpy as np
from filterpy.common import Q_continuous_white_noise


#? Defining the Symbolic Variables 
# positions 
x, y, z = sym.symbols('x, y, z')
# orientations 
yaw, pitch, roll = sym.symbols('yaw, pitch, roll')
# linear velocities 
u, v, w = sym.symbols('u, v, w')
# angular velocites 
p, q, r = sym.symbols('p, q, r')
# coefficient of drag 
Cd, Cd_f = sym.symbols('Cd, Cd_f')
# gravity 
g = sym.symbols('g')
# density
rho = sym.symbols('rho')
# input length
L = sym.symbols('L')
# Defining the Constants
m0, D, r_CP, r_CG, Ixx, Iyy, Izz, Sref, W_flap = sym.symbols('m0, D, r_CP, r_CG, Ixx, Iyy, Izz, Sref, W_flap')

#? EOMS
# moment arm 
moment_arm = r_CP - r_CG
# aerodynamic forces 
F = -((rho* (u**2) * Sref * Cd) / 2) - rho * (u**2) * Cd_f * L * W_flap

hdot = u
hddot = F/m0 - g 
omegadot = moment_arm * F / Ixx;

f_sym = sym.Matrix.vstack(sym.Matrix([[hdot],[hddot],[omegadot]]))

#? State Space 
# states
s = [x, u, p]
# inputs
i = [L]
# parameters: constants that could be subbed in later 
params = [Cd, Cd_f, Sref, rho, m0, g, Ixx, r_CG, r_CP, W_flap]
f = sym.lambdify(s + i + params, f_sym)

#? Constants to sub in 
m0 = 21.364                 # dry mass of the rocket 
D = 0.1056132               # rocket diameter 
r_CP = 219/100              # center of pressure location
r_CG = 167.67/100           # center of gravity location 
Ixx = 0.030245              # moment of inertia in x 
Iyy = 15.841                # moment of inertia in y 
Izz = Iyy                   # moment of inertia in z
Sref = np.pi * (D/2)**2    # aerodynamic reference area in m^2
W_flap = .03175             # width of the flap 
rho = 1.225                 # Constant density 
g = 9.81                    # gravitational constant 
Cd = .58                    # drag coefficient assuming constant 
Cd_f = 0.8               # drag coefficient of a flat plate 
#? Equilibrium 
x_e, u_e, p_e, L_e = 9144, 0.1, 0.01, 1

s_eq = [x_e, u_e, p_e]
i_eq = [L_e]
params_eq = [Cd, Cd_f, Sref, rho, m0, g, Ixx, r_CG, r_CP, W_flap]

f_eq = f(*s_eq,*i_eq,*params_eq)

#? State Space Calculation 
A_sym = f_sym.jacobian(s)
A_num = sym.lambdify(s + i + params, A_sym)
A = A_num(*s_eq,*i_eq,*params_eq)

B_sym = f_sym.jacobian(i)
B_num = sym.lambdify(s + i + params,B_sym)
B = B_num(*s_eq, *i_eq, *params_eq)

#? State Transition matrix 
start_time = 0
end_time = 30
total_steps = 10000
time = np.linspace(start_time,end_time,total_steps,endpoint=False)
dt = time[1] - time[0]
sigma = np.identity(3) + A * dt

#? Processing Noise Q (continuous white noise model)
Q = Q_continuous_white_noise(dim=3, dt=dt, spectral_density=1)

