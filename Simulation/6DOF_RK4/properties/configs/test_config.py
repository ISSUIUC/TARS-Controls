import math
import numpy as np

# ===== Simulator Configuration =====
# simulation output files
output_file = '../output/simulated_6dof.csv'
monte_carlo_output_folder = '../output/monte_carlo'

# Desired Apogee
des_apogee = 4572  # meters (15000 feet)

# ===== Motor Configuration =====
impulse = 9671.0  # Ns
motor_mass = 8.064  # Kg
delay = 60  # s
motor_lookup_file = '../../lookup/m2500.csv'

# ===== Rocket configuration =====
# RASAero Look Up
rasaero_lookup_file = '../lookup/RASAero.csv'

# Temporary C_d constant
C_d = 0.5

# Side Profile
C_d_s = 1.2

# length of rocket
l = 3.34

# Rocket mass without motor
rocket_dry_mass = 14.691

# radius of rocket
r_r = 0.0508 # in Meters (m)

# Center of mass of entire body
cm = np.array([3.34-2.31, 0., 0.])
cp = np.array([3.34-2.71, 0., 0.])

# Center of mass of rocket without motor
cm_rocket = np.array([3.34-1.86, 0., 0.])

# Center of mass of the motor
cm_motor = np.array([0.3755, 0., 0.])

# ===== Sensor parameters =====
# Root Mean Squared KX134 (High_G accel)
High_G_RMS = 1.9  # mg (milli g's)

# Root Mean Squared LSM6DSL (3D Gyro)
Gyro_RMS = 75  # mdps (millidegrees per second)

# Root Mean Squared MS5611 (Barometer)
Barometer_RMS = .012  # Pascal Conversion

# Rotation Error BNO08X
Bno_error = 2.5/3  # degrees of error

# ===== Other parameters =====
# flap max estension length (m)
max_ext_length = .0178

# flap max extension speed (m/s)
max_ext_spd = 0.001


# ===== Calculated quantities =====
# rocket mass with motor
rocket_total_mass = rocket_dry_mass + motor_mass

# area w/out flaps
A = math.pi*r_r**2

# side profile area
A_s = 2*r_r*l


