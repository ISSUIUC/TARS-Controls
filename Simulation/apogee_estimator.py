import numpy as np
import scipy
import matplotlib.pyplot as plt
import forces
import util.vectors as vct
import properties as prop

class Apogee: 
    
    def __init__(self, state, dt):
        self.state = state
        self.dt = dt


        pass

def step_v(time_stamp, flap_ext = 0):
    k1_v = state[2]
    
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
    return forces.get_force(np.array([pos, vel, ang_pos, ang_vel]), flap_ext, time_stamp)


    # state is in 6 by 3 matrix, access the position by accessing the state at 0
    
    