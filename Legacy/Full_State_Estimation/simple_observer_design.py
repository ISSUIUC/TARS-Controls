###############################################################################
# Illinois Space Society - SACup Avionics - Controls Team
#
# Simple observer design for full state estimation 
# Created 10/12/2021
#
# This is file derives a simplified state-space model that assumes no thrust or 
# roll moment along with constant air density throughout the flight. 
# 
# Authors:
# Ayberk Yaraneri
#
###############################################################################

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy import linalg
import control

from utils import *

# Suppress the use of scientific notation when printing small numbers
np.set_printoptions(suppress=True,precision=1,linewidth=150)
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

###############################################################################
# Known Physical Constants

dCP = 0.20          # distance between CG and CP
rho = 1.225         # air density
Sref = 0.0232       # reference area 
lref = 4.572        # reference length
m = 27.216          # mass 
Ixx, Iyy, Izz = 1.0, 1.0, 0.01    # moments of inertia (random numbers rn) 

# Aerodynamic moment coefficients.
C_l0, C_lpsi, C_malpha, C_mtheta, C_nbeta, C_nphi = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1

Cx, C_Ybeta, C_Zalpha = 0.1, 0.01, 0.01

###############################################################################
# Symoblic Variables 

# Position and derivatives in NED frame
x, y, z, xdot, ydot, zdot, xddot, yddot, zddot = \
    sym.symbols('x,y,z,xdot,ydot,zdot,xddot,yddot,zddot', real=True)

# Euler angles and derivatives in BODY frame [yaw, pitch, roll]
psi, theta, phi, phidot, thetadot, psidot, phiddot, thetaddot, psiddot = \
    sym.symbols('psi,theta,phi,phidot,thetadot,psidot, \
                 phiddot,thetaddot,psiddot', real=True)

# NED-to-BODY conversion quaternion components and derivatives                               
q1, q2, q3, q4, q1dot, q2dot, q3dot, q4dot = \
    sym.symbols('q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot', real=True)

###############################################################################
# Constructing Useful Vectors

r = sym.Matrix([x, y, z])                   # NED frame position vector
rdot = sym.Matrix([xdot, ydot, zdot])       # NED frame position derivative (velocity)
rddot = sym.Matrix([xddot, yddot, zddot])   # NED frame position 2nd derivative (acceleration)

q = sym.Matrix([q1, q2, q3, q4])                # NED-to-BODY conversion quaternion 
qdot = sym.Matrix([q1dot, q2dot, q3dot, q4dot]) # Derivative of ^^

#euler = QuatToEul(q)
#psi, theta, phi = euler[0], euler[1], euler[2]

#eulerdot = QuatToEul(qdot)
#psidot, thetadot, phidot = eulerdot[0], eulerdot[1], eulerdot[2]

omega = sym.Matrix([phidot, thetadot, psidot])           # Instantaneous angular velocity vector
#omegadot = sym.Matrix([psiddot, thetaddot, phiddot])    # Instantaneous angular acceleration vector

CP_vect = sym.Matrix([-dCP, 0, 0])     # BDY frame vector pointing from CG to CP

###############################################################################
# Constructing Forces and Monents + Axial Acceleration

# Velocity (performing NED->BDY rotation)
V_BDY = VecRotateQuat(rdot, q);
vx, vy, vz = V_BDY[0:3]

# Acceleration (performing NED->BDY rotation)
A_BDY = VecRotateQuat(rddot, q);
ax, ay, az = A_BDY[0:3]

# Velocity magnitudes for convenience
V = sym.sqrt(V_BDY.dot(V_BDY))
V2 = V_BDY.dot(V_BDY)

# Acceleration due to gravity (performing NED->BDY rotation)
accel_grav_NED = sym.Matrix([0, 0, 9.81])
accel_grav_BDY = VecRotateQuat(accel_grav_NED, q)

# Angle of attack alpha and side-slip angle beta
alpha = sym.atan2(vz, vx);
beta = sym.atan2(vy, sym.sqrt(vz**2 + vx**2))

# Aerodynamic force acting at the location of CP in BDY frame
force_aero_BDY = sym.Matrix([(Cx), \
                             (C_Ybeta*beta), \
                             (C_Zalpha*alpha)]) * (0.5*rho*V2*Sref)
fx, fy, fz = force_aero_BDY[0:3]

# Acceleration due to aero forces (performing BDY->NED rotation)
accel_aero_BDY = (1/m) * force_aero_BDY
accel_aero_NED = VecRotateQuatInv(accel_aero_BDY, q)

# Aerodynamic moment in BDY frame
moment_aero_BDY = CP_vect.cross(force_aero_BDY) + \
                  sym.Matrix([ C_l0 + C_lpsi*(psi*lref/(2*V)),               \
                               C_malpha*alpha + C_mtheta*(theta*lref/(2*V)), \
                               C_nbeta*beta + C_nphi*(phi*lref/(2*V))        ]) \
                               * 0.5*rho*V2*Sref*lref;
mx, my, mz = moment_aero_BDY[0:3]
                               
###############################################################################
# Constructing State Transtion Functions (i.e. the state derivative function f)

# State transition function of NED position and NED velocity
r_func = rdot
rdot_func = accel_grav_NED + accel_aero_NED

# State transition function of NED acceleration
rddot_func = sym.Matrix([ sym.S(0) ,
                          sym.S(0) ,
                          sym.S(0) ])

# State transition function of orientation quaterion
q_func = 0.5 * QuatMult(q, sym.Matrix([[omega], [0]]))

# State transition function of angular velocity vector (BDY frame)
omega_func = sym.Matrix([ (mx + (Iyy-Izz)*thetadot*psidot)/Ixx,
                          (my + (Izz-Ixx)*psidot*phidot)/Iyy,
                          (mz + (Ixx-Iyy)*phidot*thetadot)/Izz])

# TOTAL state transition function f:
f = sym.Matrix([ r_func     ,
                 rdot_func  ,
                 rddot_func ,
                 q_func     ,
                 omega_func ])


###############################################################################
# Constructing Measurement Functions (i.e. the state measurement function h)

h = sym.Matrix([ x       ,
                 y       ,
                 z       ,
                 xddot   ,
                 yddot   ,
                 zddot   ,
                 phidot  ,
                 thetadot,
                 psidot  ])

###############################################################################
# Linearization

states = \
    [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot,q1,q2,q3,q4,phidot,thetadot,psidot]

symbolic_variables = \
    [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot,q1,q2,q3,q4,phi,theta,psi,phidot,thetadot,psidot]

# Equilibrium states
x_e, y_e, z_e = 0, 0, 0
xdot_e, ydot_e, zdot_e = 0, 0, -100
xddot_e, yddot_e, zddot_e = 0, 0, 0

phi_e, theta_e, psi_e = 0, (sym.pi/2)-0.01, 0
phidot_e, thetadot_e, psidot_e = 0, 0, 0

q_eqb = EulToQuat(sym.Matrix([phi_e, theta_e, psi_e]))
q1_e = float(q_eqb[0]) 
q2_e = float(q_eqb[1])
q3_e = float(q_eqb[2])
q4_e = float(q_eqb[3])

A_lambda = sym.lambdify(symbolic_variables, f.jacobian(states))

C_lambda = sym.lambdify(symbolic_variables, h.jacobian(states))

A = A_lambda(x_e, y_e, z_e, xdot_e, ydot_e, zdot_e, \
             xddot_e, yddot_e, zddot_e,q1_e, q2_e, q3_e, q4_e, \
             phi_e, theta_e, psi_e, phidot_e, thetadot_e, psidot_e)

C = C_lambda(x_e, y_e, z_e, xdot_e, ydot_e, zdot_e, \
             xddot_e, yddot_e, zddot_e,q1_e, q2_e, q3_e, q4_e, \
             phi_e, theta_e, psi_e, phidot_e, thetadot_e, psidot_e)

print("### A MATRIX: \n")
print(A)

print("\n\n")

print("### C MATRIX: \n")
print(C)

###############################################################################
# Evaluating Observability 

O = control.obsv(A,C)

print("\n\n")
print("### OBSERVABILITY MATRIX: \n")
print(O)
print("\n")

obsv_rank = np.linalg.matrix_rank(O)

if obsv_rank == A.shape[0]:
    print("### OBSV MATRIX IS RANK DEFICIENT!")
    print("det(O) = " + str(obsv_rank))
else:
    print("### OBSV MATRIX IS FULL RANK!")
    print("det(O) = " + str(obsv_rank))


