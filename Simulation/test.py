import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix
from matplotlib import pyplot as plt
from statistics import mean

from helper_library import *

#? Constants(all SI units)
#* Total length of Rocket 
l = 3.0226
#* Rocket outer diameter 
D = 0.1056132
d_b = D
d_d = D
#* Body Tube Length 
L_b = 2.2352
#* Nose Cone Length 
L_n = 0.762
#* Fin thickness
T_f = 0.0029972
#* true length of the fin from inner to outer edge/ root chord 
L_m = 0.2032
#* Number of fins 
n = 3
#* fin platform area 
A_fp = 0.011532235
#* fin height 
d_f = 0.08255
#* Altitute 
z = 6096
#* Body frame velocity vector 
V_b = np.array([[100],[0],[0]])

#? Coefficient calculation 
Cd_friction = coef_v1.friction_drag(z, l, D, V_b)[0]
Re = coef_v1.friction_drag(z, l, D, V_b)[1]
Cd_body = coef_v1.body_drag(l,L_b, L_n, d_b, Cd_friction)
Cd_base = coef_v1.base_drag(d_b, d_d, Cd_body)
Cd_fin = coef_v1.fin_drag(T_f, L_m, n, A_fp, Cd_friction, d_f)
Cd = Cd_friction + Cd_body + Cd_base + Cd_fin

Rs = 200
Re_c = 51*((Rs*10**(-6)/l)**(-1.039))
print(Re_c)
