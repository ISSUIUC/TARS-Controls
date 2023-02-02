import math
import numpy as np
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6
cm = np.array([3., 0., 0.])
cp = np.array([1.06356, 0., 0.])

cm_rocket = np.array([3., 0., 0.])
cm_motor = np.array([1., 0., 0.])

# RASAero Look Up
rasaero_lookup_file = '../6DOF_RK4/LookUp/RASAero.csv'

# Temporary C_d constant
C_d = 0.5
# Side Profile
C_d_s = 1.2

# Motor Properties M2500:
impulse = 9671.0  # Ns
motor_mass = 8.064  # Kg
delay = 0  # s
motor_lookup_file = '../6DOF_RK4/LookUp/m2500.csv'

# # Motor Properties N5800:
# impulse = 20145.7 # Ns
# motor_mass = 14.826 # Kg
# delay = 0 # s
# motor_lookup_file = '../6DOF_RK4/LookUp/n5800.csv'

# # Motor Properties N2540:
# impulse = 17907 # Ns
# motor_mass = 10.700 # Kg
# delay = 0 # s
# motor_lookup_file = '../6DOF_RK4/LookUp/n2540.csv'

# rocket mass w/out motor
rocket_dry_mass = 20.2699
# rocket mass with motor
rocket_total_mass = rocket_dry_mass + motor_mass
# radius of rocket
r_r = 0.0508
# length of rocket
l = 3
# area w/out flaps
A = math.pi*r_r**2
# side profile area
A_s = 2*r_r*l
# flap max estension length (m)
max_ext_length = .0178

# Moment of Inertia


def I(total_mass): return np.diag([(1/2) * total_mass * r_r**2,
                                   (total_mass/12) * (l**2 + 3*r_r**2),
                                   (total_mass/12) * (l**2 + 3*r_r**2)])


def I_inv(total_mass): return np.diag([1/((1/2) * total_mass * r_r**2),
                                       1/((total_mass/12) * (l**2 + 3*r_r**2)),
                                       1/((total_mass/12) * (l**2 + 3*r_r**2))])


# Total Normal force
C_N_total = 9
# Total Axial force
C_A_total = 3


# simulation output file
output_file = '../6DOF_RK4/Output/simulated_6dof.csv'
