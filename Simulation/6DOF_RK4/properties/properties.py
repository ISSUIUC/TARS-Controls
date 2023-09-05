import math
import numpy as np
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6
# Center of mass of entire body
cm = np.array([3.34-2.31, 0., 0.])
cp = np.array([3.34-2.71, 0., 0.])

# RASAero Look Up
rasaero_lookup_file = '../lookup/RASAero.csv'

# Temporary C_d constant
C_d = 0.5
# Side Profile
C_d_s = 1.2

# Motor Properties M2500:
# Center of mass of rocket without motor
cm_rocket = np.array([3.34-1.86, 0., 0.])
# Center of mass of the motor
cm_motor = np.array([0.3755, 0., 0.])

impulse = 9671.0  # Ns
motor_mass = 8.064  # Kg
delay = 60  # s
motor_lookup_file = '../lookup/m2500.csv'

# # Motor Properties N5800:
# impulse = 20145.7 # Ns
# motor_mass = 14.826 # Kg
# delay = 0 # s
# motor_lookup_file = '../lookup/n5800.csv'

# # Motor Properties N2540:
# impulse = 17907 # Ns
# motor_mass = 10.700 # Kg
# delay = 0 # s
# motor_lookup_file = '../lookup/n2540.csv'

# flap max estension length (m)
max_ext_length = .0178

# flap max extension speed (m/s)
max_ext_spd = 0.001

# simulation output file
output_file = '../output/simulated_6dof.csv'
monte_carlo_output_folder = '../output/monte_carlo'

# Root Mean Squared KX134 (High_G accel)
High_G_RMS = 1.9  # mg (milli g's)

# Root Mean Squared LSM6DSL (3D Gyro)
Gyro_RMS = 75  # mdps (millidegrees per second)

# Root Mean Squared MS5611 (Barometer)
Barometer_RMS = .012  # Pascal Conversion

# Rotation Error BNO08X
Bno_error = 2.5/3  # degrees of error

# Desired Apogee
des_apogee = 4572  # meters (15000 feet)
