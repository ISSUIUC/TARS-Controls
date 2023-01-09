import numpy as np
import scipy
import matplotlib.pyplot as plt
import forces
import properties as prop

forces = forces.Forces()

def RK4(y0, dt, time_stamp, flap_ext=0) -> np.ndarray:
    '''
    Propogates State Matrix of rocket based on Runge-Kutta (RK4) Method

    Args:
        y0 (np.array): current state vector [3x6]
           [x: [pos, vel, accel, ang_pos, ang_vel, ang_accel],
            y: [pos, vel, accel, ang_pos, ang_vel, ang_accel],
            z: [pos, vel, accel, ang_pos, ang_vel, ang_accel]]
        dt (float): time step between each iteration in simulation
        time_stamp (float): current time stamp of rocket in simulation
    
    Returns:
        (np.array): state vector of rocket in x-axis [1x3]
    '''

    k1_v = (forces.get_force(np.array([y0[:,0], y0[:,1]]), flap_ext, time_stamp)[0]/prop.rocket_total_mass)[0]
    k2_v = step_v(y0[:,0], y0[:,1] + (dt/2)*k1_v, dt/2, time_stamp, flap_ext)
    k3_v = step_v(y0[:,0], y0[:,1] + (dt/2)*k2_v, dt/2, time_stamp, flap_ext)
    k4_v = step_v(y0[:,0], y0[:,1] + dt*k3_v, dt, time_stamp, flap_ext)

    v = (y0[:,1] + (1/6)*(k1_v+2*k2_v+2*k3_v+k4_v)*dt)

    k1_p = v.copy()
    k2_p = step_p(y0[:,0], y0[:,0] + (dt/2)*k1_v, dt/2)
    k3_p = step_p(y0[:,0], y0[:,0] + (dt/2)*k2_p, dt/2)
    k4_p = step_p(y0[:,0], y0[:,0] + dt*k3_p, dt)

    p = (y0[:,0] + (1/6)*(k1_p+2*k2_p+2*k3_p+k4_p)*dt)

    a = (forces.get_force(np.array([p, v]), flap_ext, time_stamp)[0]/prop.rocket_total_mass)
    return np.array([p, v, a])

def step_p(y0, y1, dt) -> np.ndarray:
    '''
    Calculates rate of change of position over given delta time for state propogation

    Args:
        y0 (np.array): current state vector [3x6]
        y1 (np.array): propogated state vector [3x6]
        dt (float): time step between iteration of RK4 (shorter than simulation dt)
    
    Returns:
        (np.array): rate of change of position (velocity) in form of state vector
    '''

    return (y1-y0)/dt # return slope (velocity)

def step_v(pos, vel, dt, time_stamp, flap_ext) -> np.ndarray:
    '''
    Calculates slope of v over given delta t for state propogation

    Args:
        pos (np.array): current posiiton state vector [3x1]
        vel (np.array): current velocity state vector [3x1]
        dt (float): time step between iteration of RK4 (shorter than simulation dt)
        flap_ext (float): current flap extention config
        time_stamp (float): current time stamp of rocket in simulation
    
    Returns:
        (np.array): rate of change of velocity (acceleration) in form of state vector
    '''

    return forces.get_force(np.array([pos, vel]), flap_ext, time_stamp)[0]/prop.rocket_total_mass # return slope (acceleration) 



