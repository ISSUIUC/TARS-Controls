from platform import mac_ver
from tracemalloc import start
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd
from filterpy.common import Q_continuous_white_noise

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation
import src.RASAero_lookup as rasaero
import src.altimeter as altimeter
import src.kalman_filter as kalman

from datetime import datetime

def accel(u, rho):
    r1 = u[0]
    v1 = u[1]
    
    # Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D)
    # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,RASaero,dic["CD"])
    Cd_total = 0.5

    F_a = -((rho* (v1**2) * Sref_a * Cd_total) / 2)
    accel_a = F_a/constants.m0 # Acceleration due to Aerodynamic Forces
    accel_f = accel_a - constants.g

    f1 = np.array([v1, accel_f])

    return f1

def rk4_step(state, dt, rho):
    # rk4 iteration 
    y1 = accel(state, rho)
    y2 = accel(state + 0.5*dt*y1, rho)
    y3 = accel(state + 0.5*dt*y2, rho)
    y4 = accel(state + dt*y3, rho)
    rk4_kp1 = state + dt*(y1 + 2*y2 + 2*y3 + y4)/6
    return rk4_kp1    

def rk4(state, dt):
    
    # Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D)
    
    while (state[1] > 0):
        
        # Define initial flap length at start of control time
        l = 0
        u = np.array([l])
        
        # A-priori 
        # kalman.priori(u)
        
        # grabbing the current states
        pos_f = state[0]
        pos_f_noise = altimeter.alt_noise(pos_f)
        vel_f = state[1]
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)
        
        # Total drag coefficient of airframe 
        
        
        # rk4 iteration 
        rk4_kp1 = rk4_step(state, dt, rho)
        
        # Append Values to the Arrays
        # dic["x"].append(float(pos_f))
        # dic["x_noise"].append(float(pos_f_noise))
        # dic["vel"].append(float(vel_f))
        # # dic["accel"].append(float(accel_f?))
        # dic["CD"].append(float(Cd_total))
        # dic["Sref"].append(float(Sref_a))
        # dic["time_sim"].append(float(t))
        
        # A-posteriori update
        # kalman.update(pos_f_noise, rk4_kp1[1], Sref_a, rho)
        
        state = rk4_kp1
                
    return state

init_cond = np.array([constants.x, constants.vx])

start_time = datetime.now()
rk4_final = rk4(init_cond, 0.03)
end_time = datetime.now()

print (end_time - start_time)
