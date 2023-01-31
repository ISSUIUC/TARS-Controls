# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np

class KalmanFilter:
    def __init__(self, dt, pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z):
        self.dt = dt
        self.x = np.zeros((9,1))
        self.Q = np.zeros((9,9))
        self.R = np.zeros((9,9))
        self.P = np.zeros((9,9))
        self.F = np.zeros((9,9))

        for i in range(3):
            self.F[3*i:3*i+2, 3*i:3*i+2] = [[1.0, dt, (dt**2) / 2],
                                    [0.0, 1.0, dt],
                                    [0.0, 0.0, 1.0]]
        self.Q = np.diag([Q_continuous_white_noise(3, dt, 0.1), Q_continuous_white_noise(3, dt, 0.1), Q_continuous_white_noise(3, dt, 0.1)])

        
        
        
        
