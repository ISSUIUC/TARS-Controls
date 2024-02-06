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
        pos_y (float): initial y position
        vel_y (float): initial y velocity
        pos_z (float): initial z position
        vel_z (float): initial z velocity
    """
    def __init__(self, dt, pos_x, pos_y, pos_z, phi, theta, psi, vel_x, vel_y, vel_z, phi_dot, theta_dot, psi_dot):
        self.dt = dt
        self.x_k = np.zeros((12,1))
        self.Q = np.zeros((12,12))
        self.R = np.diag([2., 1.9, 1.9, 1.9, 1.9, 1.9, 1.9])
        self.P_k = np.zeros((12,12), dtype=float)
        
        self.x_priori = np.zeros((12,1))
        self.P_priori = np.zeros((12,12))
        self.A = np.zeros((12,12))
        self.H = np.zeros((7,12))
        self.B = np.zeros((12,3))

        self.current_time = 0
        self.s_dt = dt

        self.x_k = np.array([pos_x, pos_y, pos_z, phi, theta, psi, vel_x, vel_y, vel_z, phi_dot, theta_dot, psi_dot]).T

        noise = Q_continuous_white_noise(2, dt, 13.)
        for i in range(6):
            self.Q[i,i] = noise[0,0]
            self.Q[i, i+6] = noise[0,1]
            self.Q[i+6, i] = noise[1,0]
            self.Q[i+6, i+6] = noise[1,1]

        # barometric altimeter 1 axis
        # accelerometer 3 axes
        self.H = np.array([[1,0,0,0,0,0,0,0,0,0,0,0],
                           [0,0,0,0,0,0,dt,0,0,0,0,0],
                           [0,0,0,0,0,0,0,dt,0,0,0,0],
                           [0,0,0,0,0,0,0,0,dt,0,0,0],
                           [0,0,0,0,0,0,0,0,0,1,0,0],
                           [0,0,0,0,0,0,0,0,0,0,1,0],
                           [0,0,0,0,0,0,0,0,0,0,0,1]], dtype=float)

        self.g = 9.81

        with open('kalman.csv', 'w') as f:
            f.write("current_time, x, y, z, r, p, y, vx, vy, vz, rdot, pdot, ydot\n")

    def priori(self, R: np.ndarray, T: np.ndarray, m: float, r:float, h:float):
        """Sets priori state and covariance
            Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
            But basically this is the prediction step
            
            Kalman Filter states are in world frame, so all data should be reported/converted to world frame
        Args:
            R (np.ndarray): rotation matrix from body to world frame
            T (np.ndarray): thrust in body frame
            m (float): mass
            rho (float): air density
            Cd (float): drag coefficient
        """
        pos_x, pos_y, pos_z = self.x_k[0:3]
        phi, theta, psi = self.x_k[3:6]
        vel_x, vel_y, vel_z = self.x_k[6:9]
        w_x, w_y, w_z = self.x_k[9:12]
        theta = -theta
        # psi, theta, phi = self.x_k[3:6]
        R = vct.body_to_world(phi, theta, psi)
        
        J_x = 1/2 * m * r**2
        J_y = 1/12 * m * h**2 + 1/4 * m * r**2
        J_z = J_y

        # State transition matrix, L means Lateral, A means angular (ur welcome)
        A_Lpp = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
        A_Lpv = R * self.dt
        A_App = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
        A_Apv = np.array([[np.cos(phi) * np.tan(theta) * self.dt, -np.sin(phi) * np.tan(theta) * self.dt, 1],
                          [np.sin(phi) * self.dt, -np.cos(phi) * self.dt, 0],
                         [np.cos(phi) / np.cos(theta) * self.dt, -np.sin(phi) / np.cos(theta) * self.dt, 0]])
        A_Lvv = np.array([[1, -w_x * self.dt, w_y * self.dt],
                          [-w_x * self.dt, -1, w_z * self.dt],
                          [-w_y * self.dt, w_z * self.dt, 1]])
        A_Avv = np.array([[1, (-w_z * J_z + w_z * J_y) / J_x * self.dt, 0],
                          [(w_z * J_z - J_x * w_z) / J_y * self.dt, -1, 0],
                          [(-w_y * J_y + J_x * w_y) / J_z * self.dt, 0, 1]])
        
        self.A = np.block([[A_Lpp, np.zeros((3,3)), A_Lpv, np.zeros((3,3))],
                        [np.zeros((3,3)), A_App, np.zeros((3,3)), A_Apv],
                        [np.zeros((3,3)), np.zeros((3,3)), A_Lvv, np.zeros((3,3))], 
                        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), A_Avv]])
        
        self.B = np.block([[np.zeros((3,3))],
                           [np.zeros((3,3))],
                           [R / m * self.dt] ,
                           [np.zeros((3,3))]])
        
        self.x_priori = self.A @ self.x_k + self.B @ T + np.array([0,0,0,0,0,0,-self.g,0,0,0,0,0])*self.dt
        self.P_priori = (self.A @ self.P_k @ self.A.T) + self.Q

    def update(self, bno_attitude, x_pos, x_accel, y_accel, z_accel, psi_vel, theta_vel, phi_vel):
        """Updates state and covariance
        
        Args:
            bno_attitude (tuple): (roll, pitch, yaw) in radians
            x_pos (float): x position
            x_accel (float): x acceleration
            y_accel (float): y acceleration
            z_accel (float): z acceleration
        """
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)
        acc = vct.body_to_world(*bno_attitude, np.array([x_accel, y_accel, z_accel]))
        ang_vel = vct.body_to_world(*bno_attitude, np.array([psi_vel, theta_vel, phi_vel]))
        y_k = np.array([x_pos, *acc, *ang_vel]).T

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt
        with open('kalman.csv', 'a') as f:
            f.write(f'{self.current_time}, {self.x_k[0]}, {self.x_k[1]}, {self.x_k[2]}, {self.x_k[6]}, {self.x_k[7]}, {self.x_k[8]}, {self.x_k[3]}, {self.x_k[4]}, {self.x_k[5]}\n')

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
        


