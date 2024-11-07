

# Error State Kalman Filter
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct
from numpy import sin, cos, tan, pi

# TODO: Figure out how to correctly implement 
class ErrorStateKalmanFilter:
    """ Kalman Filter for sensor accuracy determination
    
    Args:
        dt (float): time step
        
        pos_x (float): initial x position
        vel_x (float): initial x velocity
        accel (float): initial x acceleration
        pos_y (float): initial y position
        vel_y (float): initial y velocity
        accel (float): initial y acceleration 
        pos_z (float): initial z position
        vel_z (float): initial z velocity    
        accel_z (float): initial z acceleration 

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
    
    def _init_(self, dt, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, accel_x, accel_y, accel_z, 
               pos_x_err, pos_y_err, pos_z_err, vel_x_err, vel_y_err, vel_z_err, accel_x_err, accel_y_err, accel_z_err): 
        
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
        
        self.x_k = np.array([pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, accel_x, accel_y, accel_z]).T
        self.err_k = np.array([pos_x_err, pos_y_err, pos_z_err, vel_x_err, vel_y_err, vel_z_err, accel_x_err, accel_y_err, accel_z_err]).T
        
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
        
        
    def priori(self):
        """ 
        Sets priori state and covariance
        Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
        But basically predicts step
        Args:
        u (float): control input
        """
        self.x_k = self.x_k + self.err_k
        self.x_priori = self.F @ self.x_k
        self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q
    
    
    # TODO: Figure out if error ekf should follow similar update step to r_ekf
    #       Once implemented, add to navigation
    def update(self, bno_attitude, x_pos, x_accel, y_accel, z_accel):
        
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)
        acc = vct.body_to_world(*bno_attitude, np.array([x_accel, y_accel, z_accel])) + np.array([-9.81, 0, 0])
        y_k = np.array([x_pos, *acc]).T
        
        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt
        
        
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
        

