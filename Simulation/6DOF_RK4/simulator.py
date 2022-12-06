import numpy as np
import scipy
import matplotlib.pyplot as plt
import forces
import properties as prop

forces = forces.Forces()

def RK4(y0, dt, time_stamp, flap_ext=0):
    ''' 
    RK4 Method for Propogating Velocity

    y0 --> current state vector [3x6] [pos, vel, accel, ang_pos, ang_vel, ang_accel]

    x: [pos, vel, accel, ang_pos, ang_vel, ang_accel]
    y: [pos, vel, accel, ang_pos, ang_vel, ang_accel]
    z: [pos, vel, accel, ang_pos, ang_vel, ang_accel]

    dt --> time step
    flap_ext --> current flap extension config

    '''
    current_time = time_stamp*dt 

    k1_v = forces.get_force(y0[0], y0[1], flap_ext)/prop.mass
    k2_v = step_v(y0[0], y0[1] + (dt/2)*k1_v, dt/2, current_time, flap_ext)
    k3_v = step_v(y0[0], y0[1] + (dt/2)*k2_v, dt/2, current_time, flap_ext)
    k4_v = step_v(y0[0], y0[1] + dt*k3_v, dt, current_time, flap_ext)

    v = y0[1] + (1/6)*(k1_v+2*k2_v+2*k3_v+k4_v)*dt

    k1_p = v
    k2_p = step_p(y0[0], y0[0] + (dt/2)*k1_v, dt/2)
    k3_p = step_p(y0[0], y0[0] + (dt/2)*k2_p, dt/2)
    k4_p = step_p(y0[0], y0[0] + dt*k3_p, dt)

    p = y0[0] + (1/6)*(k1_p+2*k2_p+2*k3_p+k4_p)*dt

    a = forces.get_force(p, v, flap_ext)/prop.rocket_total_mass 

    return np.array([p, v, a])

def step_p(y0, y1, dt):
    '''
    y0 --> initial position
    y1 --> propogated position
    dt --> time step
    '''
    return (y1-y0)/dt # return slope (velocity)

def step_v(pos, y0, dt, current_time, flap_ext):
    '''
    pos --> current position
    y0 --> initial velocity
    dt --> time step
    '''
    return forces.get_force(pos, y0, flap_ext)/prop.rocket_total_mass # return slope (acceleration) 



