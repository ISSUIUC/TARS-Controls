# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct

class KalmanFilter_R:
    """Kalman Filter for 3D rotational estimation

    Args:
        dt (float): time step
        pos_x (float): initial x rotational position
        vel_x (float): initial x rotational velocity
        accel_x (float): initial x rotational acceleration
        pos_y (float): initial y rotational position
        vel_y (float): initial y rotational velocity
        accel_y (float): initial y rotational acceleration
        pos_z (float): initial z rotational position
        vel_z (float): initial z rotational velocity
        accel_z (float): initial z rotational acceleration
    """

    def __init__(self, dt, roll, pitch, yaw, w_x, w_y, w_z, a_x, a_y, a_z):
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

        self.x_k = np.array([roll, pitch, yaw, w_x, w_y, w_z, a_x, a_y, a_z]).T

        for i in range(3):
            self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
                                            [0.0, 1.0, dt],
                                            [0.0, 0.0, 1.0]]
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)


        # IMU 3 axes
        # a_x, w_x, a_y, w_y, a_z, w_z

        self.H = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
        
    
        


    # TODO: Fix priori and update step, make sure that names make sense, and add comments.
    def priori(self, R: np.ndarray, m:float, T:float, r:float, h:float, rho:float, cm:float, cp:float, Croll_a: float, Cpitch_a: float, Cyaw_a: float):
        J_x = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2
        J_z = J_y
        

        vel_x, vel_y, vel_z = self.x_k[1], self.x_k[4], self.x_k[7]
        vel_mag = np.linalg.norm([vel_x, vel_y, vel_z])

        roll_a = 0.5*rho*vel_mag*Croll_a*np.pi*(r)**2
        pitch_a = 0.5*rho*vel_mag*Cpitch_a*np.pi*(r)**2
        yaw_a = 0.5*rho*vel_mag*Cyaw_a*np.pi*(r)**2

        pos_roll = self.x_k[0]
        pos_pitch = self.x_k[1]
        pos_yaw = self.x_k[2]

        tV = np.array([T, 0, 0])
        arm = np.array([np.abs(cp - cm), 0, 0])

        tVec = tV @ R
        armVec = arm @ R

        thrustMoments = np.cross(armVec , tVec) 

        roll_p = thrustMoments[0]
        pitch_p = thrustMoments[1]
        yaw_p = thrustMoments[2]

        xdot = np.array([
            [pos_roll]
            [(roll_a + roll_p - pos_pitch*pos_yaw(J_z - J_y))/J_x]
            [1.0],
            [pos_pitch],
            [(pitch_a + pitch_p - pos_roll*pos_yaw(J_x - J_z))/J_y],
            [1.0],
            [pos_yaw],
            [(yaw_a + yaw_p - pos_roll*pos_pitch(J_y - J_x))/ J_z],
            [1.0]

            # roll_a = 0.5*rho*velocity^2*((cm - cp)/2)*rollmomentcoeff*rocketdia^2
            # pitch_a = 0.5*rho*velocity^2*((cm - cp)/2)*pitchmomentcoeff*rocketdia^2        
            # yaw_a = 0.5*rho*velocity^2*((cm - cp)/2)*yawmomentcoeff*rocketdia^2

            # For thrust calculations, mutiply Thrust by a rotation matrix in rollpitchyaw directions
            # then take the cross product of moment arm with thrust to find the forces in relative directions
            # To find moment arm distance, take scalar diff of cm and cp, and make a vector [x, 0, 0], then multiply by rotation matrix
            # finally take the cross product between the two and it will return roll pitch yaw thrust moment coefficients

        ])




        F = np.array([[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1],
                        [0, 0, 0, 0, -w_z*(-J_y + J_z)/J_x, -w_y*(-J_y + J_z)/J_x], [0, 0, 0, -w_z*(J_x - J_z)/J_y, 0, -w_x*(J_x - J_z)/J_y], [0, 0, 0, -w_y*(-J_x + J_y)/J_z, -w_x*(-J_x + J_y)/J_z, 0]])
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


