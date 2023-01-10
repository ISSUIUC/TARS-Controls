import math
import numpy as np
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6
cm = np.array([1., 0., 0.])
cp = np.array([.4, 0., 0.])

# Temporary C_d constant 
C_d = 0.5
# Side Profile
C_d_s = 1.2

# Motor Properties:
impulse = 17_907 # Ns
motor_mass = 16.281 # Kg
delay = 0 # s
# lookup_file = '../Lookup/Cesaroni_17907N2540-P_Trimmed.csv'
lookup_file = '../Lookup/cesaroni_n5800 copy.csv'

# rocket mass w/out motor
rocket_dry_mass = 23.782
# rocket mass with motor
rocket_total_mass = rocket_dry_mass + motor_mass
# radius of rocket
r_r = 0.0254
# length of rocket
l = 3
# area w/out flaps
A = math.pi*r_r**2
# side profile area
A_s = r_r*l

# simulation output file
output_file = 'simulated_6dof.csv'