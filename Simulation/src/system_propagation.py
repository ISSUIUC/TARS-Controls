from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false, interpolate
from mpl_toolkits import mplot3d
import pandas as pd
import time as timer
import random

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation
import src.RASAero_lookup as rasaero
import src.altimeter as altimeter
import src.accelerometer as accelerometer
import src.kalman_filter as kalman
import src.propellant_mass as prop
import src.interpolation as interp

# constructing the acceleration vector 
# ---> derivative of position and velocity
def accel(u, rho, Cd_total, Sref_a, thrust, rocket_mass, before_launch):
    r1 = u[0]
    v1 = u[1]
    
    # No acceleration before launch
    if (before_launch):
        return np.array([v1, 0]), 0
    
    # Calculate acceleration on the rocket
    F_a = -((rho* (v1**2) * Sref_a * Cd_total) / 2)
    accel_a = (thrust + F_a)/rocket_mass # Acceleration due to Aerodynamic Forces
    accel_f = accel_a - constants.g

    f1 = np.array([v1, accel_f])
    return f1, accel_f

def rk4_step(state, dt, rho, cd, sref, thrust, rocket_mass, before_launch):
    state = np.array([state[0], state[1]])
    #RK4: Use weighted average of 4 different slopes
    y1,a1 = accel(state, rho, cd, sref, thrust, rocket_mass, before_launch)
    y2,a2 = accel(state + 0.5*dt*y1, rho, cd, sref, thrust, rocket_mass, before_launch)
    y3,a3 = accel(state + 0.5*dt*y2, rho, cd, sref, thrust, rocket_mass, before_launch)
    y4,a4 = accel(state + dt*y3, rho, cd, sref, thrust, rocket_mass, before_launch)
    
    rk4_kp1 = state + dt*(y1 + 2*y2 + 2*y3 + y4)/6
    #! Test with weighted average of acceleration
    
    return rk4_kp1, a1

def rk4_inner(initial_state, dt, cd_file, poly_nothrust, poly_thrust, time, thrust_csv, start_mass, prop_mass_func):
    #* Returns Predicted Altitude
    
    # Initialize starting state and time
    curr_state = initial_state
    predicted_x_vals = np.array([])
    curr_mass = start_mass
    
    # Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D, 0)
    
    # Initialize polynomial fit
    poly = poly_nothrust

    if(curr_state[1] <= 0):
        return curr_state[0]
    
    # Simulate until apogee
    while (curr_state[1] >= 0):
        
        # grabbing the current states
        pos_f = curr_state[0]
        vel_f = curr_state[1]
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)
        mach = vel_f / atmosphere.speed_sound(pos_f)
        
        # Initialize variables
        before_launch = 0
        thrust = 0
                
        # Set thrust and drag values based on time
        if time < constants.thrust_start: #Before Launch
            before_launch = 1
        elif time < constants.thrust_end: #During Thrust
            curr_mass = constants.m0 - prop_mass_func(time)
            thrust = interp.thrust_interp(time, thrust_csv)
            poly = poly_thrust
        else: # After Burnout
            curr_mass = constants.mf             
            poly = poly_nothrust
        
        # Total drag coefficient of airframe 
        # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,cd_file,cd_list, 0, before_launch, before_burnout)
        Cd_total = np.poly1d(poly)(mach)
        
        # RK4 Update
        predicted_x_vals = np.append(predicted_x_vals, curr_state[0])
        next_state, accel_f = rk4_step(curr_state, dt, rho, Cd_total, Sref_a,thrust,curr_mass,before_launch)
        curr_state = next_state
        
        # Progress Time
        time += dt
    
    return max(predicted_x_vals)

def rk4_sim(initial_state, dt, cd_file, poly_nothrust, poly_thrust, desired_apogee, thrust_csv, prop_mass_func, delay, control=0):
    
    # Initialize dictionary to store values at all time steps
    sim_dict = {
    "x":[],
    "x_noise":[],
    "vel": [],
    "accel": [],
    "accel_noise": [],
    "CD": [],
    "Sref": [],
    "time_sim": [],
    "predict_alt": [],
    "flap_extension": [],
    "kalman_alt" : [],
    "kalman_vel" : [],
    "kalman_accel": []
    }
    
    start_time = int(round(timer.time()))

    # Initialize starting state and time
    curr_state = initial_state
    t = 0
    s_dt = dt #! Fix this (inner loop should be able to run at specified timestep)
    
    #mulitplier to millisec
    multiplier_to_millisecs = .001
    
    #white noise spectral dens
    spectral_density = 13.0
    
    #Error summation for integral
    e_sum = 0
    
    #Initialize Mass (with propellant)
    curr_mass = constants.m0

    mach_init = curr_state[1] / atmosphere.speed_sound(curr_state[0])
    pos_f_noise = altimeter.alt_noise(curr_state[0], mach_init)
    accel_f_noise = accelerometer.accelerometer_noise(curr_state[2])

    #* Kalman Filter Initialization
    # Initialize states (x), measurement function (H), Covariance [P], White Noise [Q], Measurement Noise Function [R]
    kalman.initialize(pos_f_noise, curr_state[1], accel_f_noise, s_dt)
    
    # Define max and min values for flap actuation
    l_max = conversion.ft_to_m(constants.max_flap_length/12) # .944 inch actuation length
    l_min = 0 # can't have negative actuation
    
    # Define initial flap length at start of control time
    u = 0
    
    # Limit flap actuation speed
    du_max = 0.001
    # Define Controller Gains
    kp, kI, kd = 0.00008, 0.0005, 0.0005
    
    # Get initial Apogee Prediction    
    # sim_dict["predict_alt"].append(38000) #!Fix
    time_elapsed = 0
    random_mulitple_step = 1
    kf_counter = 0
    predicted_apogee = 0
    # Simulate until apogee
    while (curr_state[1] >= 0 or t <= delay+1):
        # Random timestep between 10-60 ms
        # time_elapsed += s_dt
        # if time_elapsed >= (s_dt * 10):
        #     #timestep = random.randrange(6, 7)
        #     timestep = 10
        #     timestep *= multiplier_to_millisecs
            
        #     # A-priori (before current state is reached)
        #     kalman.priori(timestep, spectral_density)
        
        
        if (random_mulitple_step <= kf_counter):
            kalman.priori(random_mulitple_step * s_dt, spectral_density)

        # grabbing the current states
        pos_f = curr_state[0]
        vel_f = curr_state[1]
        accel_f = curr_state[2]
        mach = vel_f / atmosphere.speed_sound(pos_f)

        # Adding noise to initial states with sensor readings
        pos_f_noise = altimeter.alt_noise(pos_f, mach)
        accel_f_noise = accelerometer.accelerometer_noise(accel_f)
        
        # Density varies with altitude
        rho = atmosphere.density(pos_f)

        # Approximation - use the area of a circle for reference area
        Sref_a = rocket.sref_approx(constants.D, u)
        
        # Initialize variables
        before_launch = 0
        before_burnout = 0
        thrust = 0
        
        # Set thrust and drag values based on time
        if t < constants.thrust_start+delay: #Before Launch
            before_launch = 1
            before_burnout = 1
        elif t < constants.thrust_end+delay: #During Thrust
            curr_mass = constants.m0 - prop_mass_func(t-delay)
            thrust = interp.thrust_interp(t-delay, thrust_csv)
            before_burnout = 1
        else: # After Burnout
            curr_mass = constants.mf
            
            
        # Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,cd_file,sim_dict["CD"], u, before_launch, before_burnout)
        Cd_total = interp.cd_interpolation(pos_f, vel_f, 0, l_max, u, cd_file, before_launch, before_burnout)
             
        # rk4 iteration 
        next_state, accel_f = rk4_step(np.array([pos_f, vel_f]), dt, rho, Cd_total, Sref_a, thrust, curr_mass, before_launch)
        
        # Append Values to the Arrays
        sim_dict["x"].append(float(pos_f))
        sim_dict["x_noise"].append(float(pos_f_noise))
        sim_dict["vel"].append(float(vel_f))
        sim_dict["accel"].append(float(accel_f))
        sim_dict["accel_noise"].append(float(accel_f_noise))
        sim_dict["CD"].append(float(Cd_total))
        sim_dict["Sref"].append(float(Sref_a))
        sim_dict["time_sim"].append(float(t))
        sim_dict["flap_extension"].append(float(u))
        
        # Run inner RK4 to find predicted apogee
        #* (Use Kalman Filter Approximation for starting conditions)
        inner_dt = 0.1
        
        # Prediction Runtime check
        start = int(round(timer.time() * 1000))
        x_k = kalman.getStateEst()
        kalmanPos = x_k[0][0]
        kalmanVel = x_k[1][0]
        kalmanAccel = x_k[2][0]
        sim_dict["kalman_alt"].append(kalmanPos)
        sim_dict["kalman_vel"].append(kalmanVel)
        sim_dict["kalman_accel"].append(kalmanAccel)
        curr_state_est = np.array([kalmanPos, kalmanVel, kalmanAccel])
        predicted_apogee = rk4_inner(curr_state_est, inner_dt, cd_file, poly_nothrust, poly_thrust, t, thrust_csv, curr_mass, prop_mass_func)
        end = int(round(timer.time() * 1000)) - start
        
        # Calculate apogee errors 
        # Get instantaneous change in apogee error with respect to time
        apogee_error = predicted_apogee - desired_apogee
        # prev_apgee_error = sim_dict["predict_alt"][-1] - desired_apogee
        # dedt = (apogee_error - prev_apgee_error)/dt #! use dt for now, should be s_dt
        # e_sum = e_sum + (apogee_error * dt) #! use dt for now, should be s_dt
        
        # Append prediction to list
        sim_dict["predict_alt"].append(predicted_apogee)
        
        # Control Code
        if (control and not before_burnout and t > 15):
                
            prev_u = u
            # u = kp*apogee_error + kI*e_sum + kd*dedt
            u = kp*apogee_error
            # u = l_max
            u = u + np.sign((u - prev_u)/dt)*min(abs((u - prev_u)/dt), du_max)*dt
            
            # * Control Input Damping 
            if (u > l_max):
                u = l_max
            elif (u < l_min):
                u = l_min
        
        #TODO: Add check for only updating depending on s_dt
        # A-posteriori update (after current state is reached)
        '''
        if time_elapsed >= (s_dt * 10):
            kalman.update(pos_f_noise, accel_f, Sref_a, rho)
            time_elapsed = 0
        '''
        if (random_mulitple_step <= kf_counter):
            kalman.update(pos_f_noise, accel_f, Sref_a, rho)
            kf_counter = 1 
            random_mulitple_step = random.randint(1, 6)
        else:
            kf_counter += 1

        curr_state = np.array([next_state[0], next_state[1], accel_f])
        t += dt     
        
    end_time = int(round(timer.time()))
    sim_time = end_time - start_time
    return t, kalman.kalman_dict, sim_time, sim_dict
        
        
