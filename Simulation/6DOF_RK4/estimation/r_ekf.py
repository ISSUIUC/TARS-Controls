# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct

class KalmanFilter_R:
    def __init__(self, dt, pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z):
        self.dt = dt
        self.x_k = np.zeros((9,1))
        self.Q = np.zeros((9,9))
        self.R = np.diag([1.9, 1.9, 1.9, 1.9, 1.9, 1.9])
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
            # Q = np.array([[(dt**5)/20., ((dt**4)/8.)*80., (dt**3)/6.],
            #               [((dt**4)/8.)*80., ((dt**3)/8.), (dt**2)/2.],
            #               [(dt**3)/6., (dt**2)/2., dt]])
            # self.Q[3*i:3*i+3, 3*i:3*i+3] = Q*13.
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)


        # accelerometer 3 axes
        # barometric altimeter 1 axis
        self.H = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])

        self.H = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])

    def priori(self, u=0):
        self.x_priori = self.F @ self.x_k
        self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q

    def update(self, vel_x,  vel_y,  vel_z, x_accel, y_accel, z_accel):
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)
        body_rot_rate = np.array([vel_x,vel_y,vel_z])
        world_rot_rate = vct.body_to_world(self.x_k[0], self.x_k[3], self.x_k[6], body_rot_rate)
        yaw = np.arctan2(y_accel,x_accel)
        pitch = np.arctan2(z_accel, np.sqrt(y_accel**2 + x_accel**2))
        y_k = np.array([0, world_rot_rate[0], pitch, world_rot_rate[1], yaw, world_rot_rate[2]]).T

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt

    def get_state(self):
        return self.x_k

    def reset_lateral_pos(self):
        self.x_k[1:3] = [0,0]
        


