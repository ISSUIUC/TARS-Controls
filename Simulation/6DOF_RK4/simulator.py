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
        y0 (np.array): current state vector [6x3]
             x:         y:         z:
           [[pos,       pos,       pos],
            [vel,       vel,       vel],
            [accel,     accel,     accel],
            [ang_pos,   ang_pos,   ang_pos],
            [ang_vel,   ang_vel,   ang_vel],
            [ang_accel, ang_accel, ang_accel]]

        dt (float): time step between each iteration in simulation
        time_stamp (float): current time stamp of rocket in simulation
    
    Returns:
        (np.array): state vector of rocket in x-axis [6x3]
    '''

    # k1_v = forces.get_force(y0[0], y0[1], flap_ext)/prop.mass
    k1_v = (forces.get_force(np.array([y0[0], y0[1], y0[3], y0[4]]), flap_ext, time_stamp)[0]/prop.rocket_total_mass)
    k2_v = step_v(y0[0], y0[1] + (dt/2)*k1_v, y0[3], y0[4], dt/2, time_stamp, flap_ext)
    k3_v = step_v(y0[0], y0[1] + (dt/2)*k2_v, y0[3], y0[4], dt/2, time_stamp, flap_ext)
    k4_v = step_v(y0[0], y0[1] + dt*k3_v, y0[3], y0[4], dt, time_stamp, flap_ext)

    v = (y0[1] + (1/6)*(k1_v+(2*k2_v)+(2*k3_v)+k4_v)*dt)

    k1_p = y0[1].copy()
    k2_p = step_p(y0[0], y0[0] + (dt/2)*k1_p, dt/2)
    k3_p = step_p(y0[0], y0[0] + (dt/2)*k2_p, dt/2)
    k4_p = step_p(y0[0], y0[0] + dt*k3_p, dt)

    p = (y0[0] + (1/6)*(k1_p+(2*k2_p)+(2*k3_p)+k4_p)*dt)
    # a = F/m
    a = (forces.get_force(np.array([p, v, y0[3], y0[4]]), flap_ext, time_stamp)[0]/prop.rocket_total_mass)
    moment = (forces.get_force(np.array([p, v, y0[3], y0[4]]), flap_ext, time_stamp)[1]/prop.rocket_total_mass)
    alpha = moment
    # print(p)
    return np.array([p, v, a, y0[3], y0[4], alpha])

def step_p(y0, y1, dt):
    '''
    Calculates rate of change of position over given delta time for state propogation

    Args:
        y0 (np.array): current state vector [6x3]
        y1 (np.array): propogated state vector [6x3]
        dt (float): time step between iteration of RK4 (shorter than simulation dt)
    
    Returns:
        (np.array): rate of change of position (velocity) in form of state vector
    '''

    return (y1-y0)/dt # return slope (velocity)

def step_v(pos, vel, ang_pos, ang_vel, dt, time_stamp, flap_ext):
    '''
    Calculates slope of v over given delta t for state propogation

    Args:
        pos (np.array): current posiiton state vector [1x3]
        vel (np.array): current velocity state vector [1x3]
        dt (float): time step between iteration of RK4 (shorter than simulation dt)
        flap_ext (float): current flap extention config
        time_stamp (float): current time stamp of rocket in simulation
    
    Returns:
        (np.array): rate of change of velocity (acceleration) in form of state vector
    '''
    return forces.get_force(np.array([pos, vel, ang_pos, ang_vel]), flap_ext, time_stamp)[0]/prop.rocket_total_mass # return slope (acceleration) 


