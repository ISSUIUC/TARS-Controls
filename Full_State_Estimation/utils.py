###############################################################################
# Illinois Space Society - IREC 2021 Avionics Team
#
# Utility functions for Kalman filter derivation 
# Created 07/31/2021
# 
# Authors:
# Ayberk Yaraneri
#
###############################################################################

import sympy as sym
from scipy import linalg
import math

# convert distance from feet to meters
def ft_to_m(measurement):
    return (measurement / 3.2808) 

# Optimize a linear quadratic regulator by solving the ricatti equation
def lqr(A, B, Q, R):
    P = linalg.solve_continuous_are(A, B, Q, R)
    K = linalg.inv(R) @  B.T @ P
    return K

# Perform quaternion mutliplication
def QuatMult(qA, qB):

    mat = sym.Matrix([[ qA[3],-qA[2], qA[1], qA[0] ],
                      [ qA[2], qA[3],-qA[0], qA[1] ],
                      [-qA[1], qA[0], qA[3], qA[2] ],
                      [-qA[0],-qA[1],-qA[2], qA[3] ]])

    return mat * qB;

# Rotate vector v using quaternion q
def VecRotateQuat(v, q):

    qInv = sym.Matrix([-q[0], -q[1], -q[2], q[3]])

    quatVect = sym.Matrix([[v], [0]])

    quatOut = QuatMult(QuatMult(q, quatVect), qInv)

    vecOut = sym.Matrix(quatOut[0:3])

    return vecOut

# Rotate vector v using inverse of quaternion q
def VecRotateQuatInv(v, q):

    qInv = sym.Matrix([-q[0], -q[1], -q[2], q[3]])

    quatVect = sym.Matrix([[v], [0]])

    quatOut = QuatMult(QuatMult(qInv, quatVect), q)

    vecOut = sym.Matrix(quatOut[0:3])

    return vecOut

# Convert from a 321 Euler rotation sequence in radians to a Quaternion
def EulToQuat(Euler):

    Euler    = Euler * 0.5;
    cosPhi   = sym.cos(Euler[0]);
    sinPhi   = sym.sin(Euler[0]);
    cosTheta = sym.cos(Euler[1]);
    sinTheta = sym.sin(Euler[1]);
    cosPsi   = sym.cos(Euler[2]);
    sinPsi   = sym.sin(Euler[2]);

    quatOut = sym.Matrix([[ sinPhi*cosTheta*cosPsi - cosPhi*sinTheta*sinPsi ],
                          [ cosPhi*sinTheta*cosPsi + sinPhi*cosTheta*sinPsi ],
                          [ cosPhi*cosTheta*sinPsi - sinPhi*sinTheta*cosPsi ],
                          [ cosPhi*cosTheta*cosPsi + sinPhi*sinTheta*sinPsi ]])

    return quatOut

# TODO: Fix quaternion convention used in this function
# DO NOT USE - Uses wrong quaternion convention of w,x,y,z. Should be x,y,z,w
# Convert from a quaternion to a 321 Euler rotation sequence in radians
#def QuatToEul(q):
#
#    q0,q1,q2,q3 = q[0],q[1],q[2],q[3]
#
#    phi = sym.atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1**2 + q2**2))
#
#    theta = sym.asin(2*(q0*q2 - q3*q1))
#
#    psi = sym.atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
#
#    return sym.Matrix([[psi]   ,
#                       [theta] ,
#                       [phi]   ])
