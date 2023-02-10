# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np

class KalmanFilter:
    def __init__(self, dt, pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z):
        self.dt = dt
        self.x_k = np.zeros((9,1))
        self.Q = np.zeros((9,9))
        self.R = np.diag([1., 100., 100., 100.])
        self.P_k = np.zeros((9,9))
        self.x_priori = np.zeros((9,1))
        self.P_priori = np.zeros((9,9))
        self.F = np.zeros((9,9))
        self.H = np.zeros((4,9))

        self.current_time = 0
        self.s_dt = dt

        self.x_k = np.array([pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z]).T

        for i in range(3):
            self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
                                    [0.0, 1.0, dt],
                                    [0.0, 0.0, 1.0]]
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 0.1)

        # accelerometer 3 axes
        # barometric altimeter 1 axis
        self.H = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

    def priori(self, u):
        self.x_priori = self.F @ self.x_k
        self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q

    def update(self, x_pos, x_accel, y_accel, z_accel):
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)
        
        y_k = np.array([x_pos, x_accel, y_accel, z_accel]).T

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt

    def get_state(self):
        return self.x_k
        
        
        
