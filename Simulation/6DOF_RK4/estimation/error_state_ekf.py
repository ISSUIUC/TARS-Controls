# error state kalman fiter
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct


# TODO: Figure out how to correctly implement 
class ErrorStateKalmanFilter:
    """ Kalman Filter for sensor accuracy determination
    
    Args:
        dt (float): time step
        pos_x_err (float): initial x position error
        vel_x_err (float): initial x velocity error
        accel_x_err (float): initial x acceleration error
        pos_y_err (float): initial y position error
        vel_y_err (float): initial y velocity error
        accel_y_err (float): initial y acceleration error
        pos_z_err (float): initial z position error
        vel_z_err (float): initial z velocity error     
        accel_z_err (float): initial z acceleration error

    """
    
    def _init_(self, dt, pos_x_err, pos_y_err, pos_z_err, vel_x_err, vel_y_err, vel_z_err, accel_x_err, accel_y_err, accel_z_err): 
        
        self.dt = dt
        
        self.x_k = np.zeros((9,1))
        self.Q = np.zeros((9,9))
        self.R = np.diag([2., 1.9, 1.9, 1.9])
        self.P_k = np.zeros((9,9))
        self.x_priori = np.zeros((9,1)) 
        self.P_priori = np.zeros((9,9))
        self.F = np.zeros((9,9))
        self.H = np.zeros((4,9))
    
        self.current_time = 0 
        self.s_dt = dt
        
        self.x_k = np.array([pos_x_err, pos_y_err, pos_z_err, vel_x_err, vel_y_err, vel_z_err, accel_x_err, accel_y_err, accel_z_err]).T
        
        # Initialize F
        for i in range(3):
            self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
                                            [0.0, 1.0, dt],
                                            [0.0, 0.0, 1.0]]
        # Initialize Q and R covariance matrices 
        self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)
        
        self.H = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
        
        
        def priori():
            """ 
            Sets priori state and covariance
            Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
            But basically predicts step
            Args:
            u (float): control input
            """
            # TODO: figure out if self.F has to be jacobian(f w/ respect to x) for ekf: f(x_{k-1})
            self.x_priori = self.F @ self.x_k
            self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q
        
        
    def get_state(self):
        """Returns current state
        
        Returns:
            np.array: current state
        """
        return self.x_k
    
    def get_covariance(self):
        """Returns current covariance
        
        Returns:
            np.array: current covariance
        """
        return np.diag(self.P_k)
    
    def reset_lateral_pos(self):
        """Resets lateral position to 0
        """
        self.x_k[1:3] = [0,0]
