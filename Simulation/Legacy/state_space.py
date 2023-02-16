import time
import matplotlib.pyplot as plt
from scipy import linalg
import sympy as sym
import numpy as np
import src.rotation as rotation

#? Defining the Constants
m0, D, r_CP, r_CG, Ixx, Iyy, Izz, Sref, W_flap = sym.symbols('m0, D, r_CP, r_CG, Ixx, Iyy, Izz, Sref, W_flap')

def rotate_yaw(phe):
    # rotation about the z-axis, takes in angle phe in the function (rad)
    return sym.Matrix([[sym.cos(phe), -sym.sin(phe), 0],[sym.sin(phe), sym.cos(phe), 0],[0, 0, 1]])

def rotate_pitch(theta):
    # rotation about the y-axis, takes in angle theta in the function (rad)
    return sym.Matrix([[sym.cos(theta), 0, sym.sin(theta)],[0, 1, 0],[-sym.sin(theta), 0, sym.cos(theta)]])

def rotate_roll(phi):
    # rotation about the x-axis, take sin angle phi in the function (rad)
    return sym.Matrix([[1, 0, 0],[0, sym.cos(phi), -sym.sin(phi)], [0, sym.sin(phi), sym.cos(phi)]])

def body_aero(velocity_body):
    #TODO: Disperse comments to properly explain the lines of code (atmosphere.py for reference)
    # Takes in np.array velocity_body 
    # rotation matrix that transforms the body frame to the aerodyanmic frame 
    # takes in beta (side-slip angle), and alpha (angles of attack)
    # needs the velocities in the body frame to calculate beta and alpha
    vx_b = velocity_body[0]
    vy_b = velocity_body[1]
    vz_b = velocity_body[2]

    alpha = sym.atan2(vz_b,vx_b)
    beta = sym.atan2(vy_b, sym.sqrt(vx_b**2 + vz_b**2))
    R_ba = sym.Matrix([[sym.cos(beta)*sym.cos(alpha), sym.sin(beta), sym.cos(beta)*sym.sin(alpha)], [-sym.sin(beta)*sym.cos(alpha), sym.cos(beta), 
    -sym.sin(beta)*sym.sin(alpha)],[-sym.sin(alpha), 0, sym.cos(alpha)]])
    return R_ba

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
Cd = sym.symbols('Cd')
# gravity 
g = sym.symbols('g')
# density
rho = sym.symbols('rho')
# input length
L = sym.symbols('L')

#? turning variables into matrices 
pos_f = sym.Matrix([[x], [y], [z]])
or_f = sym.Matrix([[yaw], [pitch], [roll]])
vel_f = sym.Matrix([[u], [v], [w]])
ang_vel_f = sym.Matrix([[p], [q], [r]])
I_b = sym.Matrix([[Ixx, 0, 0], [0, Iyy, 0], [0, 0, Izz]])

#Symbolic Rotation Matrices
R_fb = rotate_yaw(yaw) @ rotate_pitch(pitch) @ rotate_roll(roll)
R_bf = R_fb.inv()

# Fixed frame to body frame
vel_b = R_bf @ vel_f
# Body to Aerodynamic
R_ba = body_aero(vel_b)
R_ab = R_ba.T #rotation matricies are orthonormal -> transpose  = inverse but takes less computation time
vel_a = R_ab @ vel_b

# Moment arm
moment_arm_b = sym.Matrix([[r_CP-r_CG],[0],[0]])
moment_arm_a = R_ab @ moment_arm_b

#  Magnitude of velocity
v_mag = vel_a.norm()

# Aerodynamic forces
F_a = -((rho* (v_mag**2) * Sref * Cd) / 2) * (vel_a/v_mag)

# Linear Acceleration in the aerodynamic frame
accel_a = F_a/m0

# Torque
torque_a = moment_arm_a.cross(F_a)
torque_b = R_ba @ torque_a

# Angular acceleration in the body frame
ang_accel_b = I_b.inv() @ torque_b

# Angular acceleration in the fixed frame
ang_accel_f = R_fb @ ang_accel_b

# Linear Acceleration in the body frame
accel_b = R_ba @ accel_a

# Linear Acceleration in the fixed frame
accel_f = R_fb @ accel_b - sym.Matrix([[g],[0],[0]])

# EOMs
hdot = vel_f
hddot = accel_f
or_dot = ang_vel_f
or_ddot = ang_accel_f
f_sym = sym.Matrix.vstack(sym.Matrix([[hdot],[hddot],[or_dot],[or_ddot]]))

# State Space
s = [x,y,z,u,v,w,yaw,pitch,roll,p,q,r]

i = [L]

params = [Cd, Sref, rho, W_flap, m0, g, Ixx, Iyy, Izz, D, r_CG, r_CP]

f = sym.lambdify(s+i+params, f_sym)

m0 = 21.364                 # dry mass of the rocket 
D = 0.1056132               # rocket diameter 
r_CP = 219/100              # center of pressure location
r_CG = 167.67/100           # center of gravity location 
Ixx = 0.030245              # moment of inertia in x 
Iyy = 15.841                # moment of inertia in y 
Izz = Iyy                   # moment of inertia in z
Sref = sym.pi * (D/2)**2    # aerodynamic reference area in m^2
W_flap = .03175
rho = 1.225
g = 9.81
Cd = .58

u_e, v_e, w_e = .1,.1,.1
x_e, y_e, z_e = 9144,0,0
yaw_e, pitch_e, roll_e = 0.1,0.1,0.1
p_e, q_e, r_e = 0.01,0.01,0.01
L_e = 1

s_eq = [x_e, y_e, z_e, u_e, v_e, w_e, yaw_e, pitch_e, roll_e, p_e, q_e, r_e]
i_eq = [L_e]
params_eq = [Cd, Sref, rho, W_flap, m0, g, Ixx, Iyy, Izz, D, r_CG, r_CP]

f_eq = f(*s_eq,*i_eq,*params_eq)

A_sym = f_sym.jacobian(s)
A_num = sym.lambdify(s+i+params, A_sym)
A = A_num(*s_eq,*i_eq,*params_eq)
print(A)