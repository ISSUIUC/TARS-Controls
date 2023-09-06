import math
import numpy as np

# ===== Motor Configuration =====
impulse = 9671.0  # Ns
motor_mass = 8.064  # Kg
delay = 60  # s
motor_lookup_file = '../../lookup/m2500.csv'

# ===== Rocket configuration =====
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

# ===== Other parameters =====
# flap max estension length (m)
max_ext_length = .0178


# ===== Calculated quantities =====
# rocket mass with motor
rocket_total_mass = rocket_dry_mass + motor_mass

# area w/out flaps
A = math.pi*r_r**2

# side profile area
A_s = 2*r_r*l


