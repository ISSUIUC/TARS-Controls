

# all the linear positions
from filterpy.common import Q_continuous_white_noise
import numpy as np
import util.vectors as vct
from numpy import sin
from numpy import cos
from numpy import tan
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
        self.H = np.zeros((7,12))

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
        # gyro 3 axes
        self.H = np.array([[1,0,0,0,0,0,0,0,0,0,0,0],
                           [0,0,0,0,0,0,0,0,dt,0,0,0],
                           [0,0,0,0,0,0,0,dt,0,0,0,0],
                           [0,0,0,0,0,0,dt,0,0,0,0,0],
                           [0,0,0,0,0,0,0,0,0,0,0,1],
                           [0,0,0,0,0,0,0,0,0,0,1,0],
                           [0,0,0,0,0,0,0,0,0,1,0,0]], dtype=float)

        self.g = 9.81

        with open('kalman.csv', 'w') as f:
            f.write("current_time,x,y,z,r,p,y,vx,vy,vz,rdot,pdot,ydot\n")

    def priori(self, R: np.ndarray, T: float, m: float, r:float, h:float):
        """Sets priori state and covariance
            Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
            But basically this is the prediction step
            
            Kalman Filter states are in world frame, so all data should be reported in/converted to world frame
        Args:
            R (np.ndarray): rotation matrix from body to world frame
            T (float): thrust magnitude (assumed to be along rocket body)
            m (float): mass
            rho (float): air density
            Cd (float): drag coefficient
        """
        # Transformation array switching x and z and the negative of y
        # P converts from the state in Anshuk's paper (z vertical up) to the state in our body frame (x vertical up) to the world frame
        
        P = np.array([[0, 0, 1],
                      [0, -1, 0],
                      [1, 0, 0]])
        transformation = np.block([[R@P, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
                                   [np.zeros((3,3)), R@P, np.zeros((3,3)), np.zeros((3,3))],
                                   [np.zeros((3,3)), np.zeros((3,3)), R@P, np.zeros((3,3))],
                                   [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), R@P]])
        
        x_k = np.linalg.inv(transformation) @ self.x_k  
        pos_x, pos_y, pos_z = x_k[0:3]
        phi, theta, psi = x_k[3:6]
        vel_x, vel_y, vel_z = x_k[6:9]
        w_x, w_y, w_z = x_k[9:12]
        # R = vct.body_to_world(phi, theta, psi)
        
        J_z = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2
        J_x = J_y
        # J_z = J_y
        # A_Lpp = np.eye(3)
        # A_Lpv = np.eye(3) @ R * self.dt
        # A_App = np.eye(3)
        # A_Apv = np.array([[0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)],
        #                   [0, np.cos(phi), -np.sin(phi)],
        #                   [1 / self.dt, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)]]) * self.dt
        # A_Lvv = np.array([[1, w_z * self.dt, -w_y * self.dt],
        #                   [-w_z * self.dt, 1, w_x * self.dt],
        #                   [w_y * self.dt, -w_x * self.dt, 1]])
        # A_Avv = np.array([[1, w_z * J_y / J_x * self.dt, -w_y * J_z / J_x * self.dt],
        #                   [-w_z * J_x / J_y * self.dt, 1, w_x * J_z / J_y * self.dt],
        #                   [w_y * J_x / J_z * self.dt, w_x * J_y / J_z * self.dt, 1]])
        
        # A = np.block([[A_Lpp, np.zeros((3,3)), A_Lpv, np.zeros((3,3))],
        #               [np.zeros((3,3)), A_App, np.zeros((3,3)), A_Apv],
        #               [np.zeros((3,3)), np.zeros((3,3)), A_Lvv, np.zeros((3,3))], 
        #               [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), A_Avv]])
        g = 9.81
        x_dot = np.array([
            vel_x*cos(psi)*cos(theta) + vel_y*(sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi)) + vel_z*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)),
            vel_x*sin(psi)*cos(theta) + vel_y*(sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi)) + vel_z*(-sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)),
            -vel_x*sin(theta) + vel_y*sin(phi)*cos(theta) + vel_z*(cos(phi)*cos(theta)),
            (w_y*sin(phi))/(cos(theta)) + (w_z*cos(phi))/(cos(theta)),
            w_y*cos(phi) - w_z*sin(phi),
            w_x + w_y*sin(phi)*tan(theta) + w_z*cos(phi)*tan(theta),
            (g*m*sin(theta)+m*vel_y*w_z-m*vel_z*w_y) / m,
            (-g*m*sin(phi)*cos(theta) - m*vel_x*w_z+m*vel_z*w_x) / m,
            (abs(T) - g*m*cos(phi)*cos(theta) + m*vel_x*w_y - m*vel_y*w_x) / m,
            (J_y*w_z*w_y - J_z*w_z*w_y) / J_x,
            (-J_x*w_z*w_x + J_z*w_z*w_x) / J_y,
            (J_x*w_y*w_x - J_y*w_x*w_y) / J_z
        ])
        # B = np.block([[np.zeros((3,3))],
        #               [np.zeros((3,3))],
        #               [np.eye(3) / m * self.dt],
        #               [np.zeros((3,3))]])
        # Puts A and B matrix back into our states
        # For transformation => x = transformation inverse (xbar) transformation (A matrix)
        # (B matrix) = transformation inverse (b bar)
        # change of state variables on page 50-->ECE 515 Course Notes https://arxiv.org/abs/2007.01367
        # A = np.linalg.inv(transformation) @ A @ transformation
        # B = np.linalg.inv(transformation) @ B
        # g = np.array([0,0,0,0,0,0,self.g*np.sin(theta), -self.g*np.sin(theta)*np.cos(theta), -self.g*np.cos(theta)*np.cos(phi),0,0,0])*self.dt
        # TODO: Check if transformations are correct (should there be an inverse behind x_dot) jk it seems right now
        self.x_priori = self.x_k + (transformation@x_dot)*self.dt
        
        A = transformation@np.array([[0, 0, 0, vel_y*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)) + vel_z*(-sin(phi)*sin(theta)*cos(psi) + sin(psi)*cos(phi)), -vel_x*sin(theta)*cos(psi) + vel_y*sin(phi)*cos(psi)*cos(theta) + 
vel_z*cos(phi)*cos(psi)*cos(theta), -vel_x*sin(psi)*cos(theta) + vel_y*(-sin(phi)*sin(psi)*sin(theta) - cos(phi)*cos(psi)) + vel_z*(sin(phi)*cos(psi) - sin(psi)*sin(theta)*cos(phi)), cos(psi)*cos(theta), sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi), sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi), 0, 0, 0], [0, 0, 0, vel_y*(-sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)) + vel_z*(-sin(phi)*sin(psi)*sin(theta) - cos(phi)*cos(psi)), -vel_x*sin(psi)*sin(theta) + vel_y*sin(phi)*sin(psi)*cos(theta) + vel_z*sin(psi)*cos(phi)*cos(theta), vel_x*cos(psi)*cos(theta) + vel_y*(sin(phi)*sin(theta)*cos(psi) - sin(psi)*cos(phi)) + vel_z*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)), sin(psi)*cos(theta), sin(phi)*sin(psi)*sin(theta) + cos(phi)*cos(psi), -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi), 0, 0, 0], [0, 0, 0, vel_y*cos(phi)*cos(theta) - vel_z*sin(phi)*cos(theta), -vel_x*cos(theta) - vel_y*sin(phi)*sin(theta) - vel_z*sin(theta)*cos(phi), 0, -sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta), 0, 0, 0], [0, 0, 0, w_y*cos(phi)/cos(theta) - w_z*sin(phi)/cos(theta), w_y*sin(phi)*sin(theta)/cos(theta)**2 + w_z*sin(theta)*cos(phi)/cos(theta)**2, 0, 0, 0, 0, 0, sin(phi)/cos(theta), cos(phi)/cos(theta)], [0, 0, 0, -w_y*sin(phi) - w_z*cos(phi), 0, 0, 0, 0, 0, 0, cos(phi), -sin(phi)], [0, 0, 0, w_y*cos(phi)*tan(theta) - w_z*sin(phi)*tan(theta), w_y*(tan(theta)**2 + 1)*sin(phi) + w_z*(tan(theta)**2 + 1)*cos(phi), 0, 0, 0, 0, 1, sin(phi)*tan(theta), cos(phi)*tan(theta)], [0, 0, 0, 0, g*cos(theta), 0, 0, w_z, -w_y, 0, -vel_z, vel_y], [0, 0, 0, -g*cos(phi)*cos(theta), g*sin(phi)*sin(theta), 
0, -w_z, 0, w_x, vel_z, 0, -vel_x], [0, 0, 0, g*sin(phi)*cos(theta), g*sin(theta)*cos(phi), 0, w_y, -w_x, 0, -vel_y, vel_x, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (J_y*w_z - J_z*w_z)/J_x, (J_y*w_y - J_z*w_y)/J_x], [0, 0, 0, 0, 0, 0, 0, 0, 0, (-J_x*w_z + J_z*w_z)/J_y, 0, (-J_x*w_x + J_z*w_x)/J_y], [0, 0, 0, 0, 0, 0, 0, 0, 0, (J_x*w_y - J_y*w_y)/J_z, (J_x*w_x - J_y*w_x)/J_z, 0]])
        
        self.P_priori = (A @ self.P_k @ A.T) + self.Q

    def update(self, bno_attitude, x_pos, x_accel, y_accel, z_accel, phi_vel, theta_vel, psi_vel):
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
        ang_vel = vct.body_to_world(*bno_attitude, np.array([phi_vel, theta_vel, psi_vel]))
        y_k = np.array([x_pos, *acc, *ang_vel]).T

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt
        with open('kalman.csv', 'a') as f:
            f.write(f'{self.current_time}, {self.x_k[0]}, {self.x_k[1]}, {self.x_k[2]}, {self.x_k[3]}, {self.x_k[4]}, {self.x_k[5]}, {self.x_k[6]}, {self.x_k[7]}, {self.x_k[8]}\n')

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
        


