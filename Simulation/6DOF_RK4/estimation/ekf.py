# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct

class KalmanFilter:
    """Kalman Filter for 3D position estimation
    
    Args:
        dt (float): time step
        pos_x (float): initial x position
        vel_x (float): initial x velocity
        accel_x (float): initial x acceleration
        pos_y (float): initial y position
        vel_y (float): initial y velocity
        accel_y (float): initial y acceleration
        pos_z (float): initial z position
        vel_z (float): initial z velocity
        accel_z (float): initial z acceleration
    """
    def __init__(self, dt, pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z):
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

        for i in range(3):
            # self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
            #                                 [0.0, 1.0, dt],
            #                                 [0.0, 0.0, 1.0]]
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)


        # barometric altimeter 1 axis
        # accelerometer 3 axes
        self.H = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

    def priori(self, R: np.ndarray, T: np.ndarray, m: float, rho: float, Cd: float):
        """Sets priori state and covariance
            Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
            But basically this is the prediction step
        Args:
            R (np.ndarray): rotation matrix from body to world frame
            T (np.ndarray): thrust in body frame
            m (float): mass
            rho (float): air density
            Cd (float): drag coefficient
        """
        # State transition matrix
        A = np.block([[np.eye(3), self.dt*np.eye(3), (self.dt**2)/2*np.eye(3)],
                      [np.zeros((3,3)), np.eye(3), self.dt*np.eye(3)],
                      [np.zeros((3,3)), np.zeros((3,3)), np.eye(3)]])
        
        # B matrix (generally used for control input but we're using thrust as "control" input)
        B = np.block([[np.zeros((3,3))],
                      [np.zeros((3,3))],
                      [R/m]])
        
        self.x_priori = A @ self.x_k + B @ T + np.array([0,0,0,0,0,0,-9.81,0,0])
        self.P_priori = (A @ self.P_k @ A.T) + self.Q

    def update(self, bno_attitude, x_pos, x_accel, y_accel, z_accel):
        """Updates state and covariance
        
        Args:
            bno_attitude (tuple): (roll, pitch, yaw) in radians
            x_pos (float): x position
            x_accel (float): x acceleration
            y_accel (float): y acceleration
            z_accel (float): z acceleration
        """
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
        


