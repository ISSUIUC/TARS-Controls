from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd
import time as timer

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

# constructing the acceleration vector 
# ---> derivative of position and velocity
def accel(u, rho, Cd_total, Sref_a):
    r1 = u[0]
    v1 = u[1]
    
    F_a = -((rho* (v1**2) * Sref_a * Cd_total) / 2)
    accel_a = F_a/constants.m0 # Acceleration due to Aerodynamic Forces
    accel_f = accel_a - constants.g

    f1 = np.array([v1, accel_f])

    return f1

def rk4_step(state, dt, rho, cd, sref):
    # rk4 iteration 
    y1 = accel(state, rho, cd, sref)
    y2 = accel(state + 0.5*dt*y1, rho, cd, sref)
    y3 = accel(state + 0.5*dt*y2, rho, cd, sref)
    y4 = accel(state + dt*y3, rho, cd, sref)
    rk4_kp1 = state + dt*(y1 + 2*y2 + 2*y3 + y4)/6
    return rk4_kp1

def rk4_inner(initial_state, dt, cd_file, poly):
    #* Returns Predicted Altitude
    
    # Initialize starting state and time
    curr_state = initial_state
    t = 0    
    predicted_x_vals = np.array([])
    
    # Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D, 0)
    
    # Simulate until apogee
    while (curr_state[1] > 0):
        
        # grabbing the current states
        pos_f = curr_state[0]
        pos_f_noise = altimeter.alt_noise(pos_f)
        vel_f = curr_state[1]
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)
        mach = curr_state[1] / atmosphere.speed_sound(curr_state[0])
        
        # Total drag coefficient of airframe 
        # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,cd_file,dict["CD"])
        # Cd_total = 0.5
        Cd_total = np.poly1d(poly)(mach)
        
        # rk4 iteration
        predicted_x_vals = np.append(predicted_x_vals, curr_state[0])
        next_state = rk4_step(curr_state, dt, rho, Cd_total, Sref_a)
        curr_state = next_state
        t += dt    
    
    return max(predicted_x_vals)

def rk4_sim(initial_state, pos_f_noise, dt, cd_file, poly, control=0):
    
    # Initialize dictionary to store values at all time steps
    sim_dict = {
    "x":[],
    "x_noise":[],
    "vel": [],
    "accel": [],
    "CD": [],
    "Sref": [],
    "time_sim": [],
    "predict_alt": [],
    "predict_update_alt": [],
    "flap_extension": []
    }
    
    start_time = int(round(timer.time()))

    # Initialize starting state and time
    curr_state = initial_state
    t = 0
    s_dt = dt #!

    #* Kalman Filter Initialization
    # Initialize states (x), measurement function (H), Covariance [P], White Noise [Q], Measurement Noise Function [R]
    kalman.initialize(pos_f_noise, curr_state[1], s_dt)
    
    # Define max and min values for flap actuation
    l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
    l_min = 0 # can't have negative actuation
    
    # Define initial flap length at start of control time
    l = 0
    u = np.array([l])
    
    # Simulate until apogee
    while (curr_state[1] > 0):
        
        # A-priori (before current state is reached)
        kalman.priori(u)
        
        # grabbing the current states
        pos_f = curr_state[0]
        pos_f_noise = altimeter.alt_noise(pos_f)
        vel_f = curr_state[1]
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)
        
        # Total drag coefficient of airframe 
        Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,cd_file,sim_dict["CD"], u[0])
        # Cd_total = 0.5

        # Approximation - use the area of a circle for reference area
        Sref_a = rocket.sref_approx(constants.D, u[0])

        # rk4 iteration 
        next_state = rk4_step(curr_state, dt, rho, Cd_total, Sref_a)
        
        # Append Values to the Arrays
        sim_dict["x"].append(float(pos_f))
        sim_dict["x_noise"].append(float(pos_f_noise))
        sim_dict["vel"].append(float(vel_f))
        # dic["accel"].append(float(accel_f?))
        sim_dict["CD"].append(float(Cd_total))
        sim_dict["Sref"].append(float(Sref_a))
        sim_dict["time_sim"].append(float(t))
        sim_dict["flap_extension"].append(float(u[0]))
        
        # Run inner RK4 to find predicted apogee
        #* (Use Kalman Filter Approximation for starting conditions)
        inner_dt = 0.1
        
        # Prediction Runtime check
        start = int(round(timer.time() * 1000))
        predicted_apogee = rk4_inner(curr_state, inner_dt, cd_file, poly)
        end = int(round(timer.time() * 1000)) - start
        sim_dict["predict_alt"].append(predicted_apogee)
        
        # Control Code
        if (control):
            k = -0.003
            u = -k*predicted_apogee
            # u = u_t1 + np.sign((u_t1 - u[0])/dt)*min(abs((u_t1 - u[0])/dt), du_max)*dt
            
            #* Control Input Damping 
            if (u > l_max):
                u = l_max
            elif (u < l_min):
                u = l_min
                
            u = [u]
            
            # print("Flap Extension: ", u[0])
            
        #TODO: Add check for only updating depending on s_dt
        # A-posteriori update (after current state is reached)
        kalman.update(pos_f_noise, curr_state[1], Sref_a, rho)
        curr_state = next_state
        t += dt     
        
    end_time = int(round(timer.time()))
    sim_time = end_time - start_time

   
    return t, kalman.kalman_dic, sim_time, sim_dict
        
        