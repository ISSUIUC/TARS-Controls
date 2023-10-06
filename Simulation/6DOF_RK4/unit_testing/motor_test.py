import numpy as np
import sys
import os
import math
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
import dynamics.forces as forces 
import dynamics.rocket as rocket_model
import environment.atmosphere as atmosphere
import dynamics.motor as mot
import util.vectors as vct
import pandas as pd
import properties.properties as prop

# Defining necessary variables
start_time = prop.delay
total_impulse = 9671.0
total_mass = 8.064
current_mass = total_mass
#git pucoast_time = delay
alignment = np.array([0, 0])
start_time = prop.delay
cur_line = 0
cm = np.array([3.34-2.31, 0., 0.])
cm_rocket = np.array([3.34-1.86, 0., 0.])
cm_motor = np.array([0.3755, 0., 0.])
rocket_dry_mass = 14.691
rocket_motor_mass = 8.064
rocket_total_mass = rocket_dry_mass + rocket_motor_mass


