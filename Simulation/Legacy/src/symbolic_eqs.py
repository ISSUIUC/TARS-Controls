from sympy import *
import numpy as np
import rotation as rotation
import constants as constants
import rocket as rocket


rho, va, vbx, vby, vbz, s_ref, cd, c_m, c_p, g, m, Ix, Iy, Iz = symbols('rho, va, vbx, vby, vbz, s_ref, cd, c_m, c_p, g, m, Ix, Iy, Iz')

vb = Matrix([[vbx], [vby], [vbz]])
or_f = constants.init_or_f
R_fb = rotation.yaw(or_f[0][0]) @ rotation.pitch(or_f[1][0]) @ rotation.roll(or_f[2][0])
R_bf = np.linalg.inv(R_fb)
R_ab = rotation.body_aero(vb)

Fd = Matrix([(rho*va*s_ref*cd)/2, 0, 0]).T
Fg = Matrix([g*m, 0, 0])
F = (R_bf @ R_ab @ Fd.T) - Fg
a_f = F/m
I = Matrix([[Ix, 0, 0],[0, Iy, 0], [0, 0, Iz]])
l_M = Matrix([c_m - c_p, 0, 0])
M_a = ((np.transpose(R_ab))@l_M).cross(Fd)
ang_accel = R_bf@R_ab@(I.inv()@M_a)
print(ang_accel)
