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
        self.P_k = np.eye(9)
        self.x_priori = np.zeros((9,1))
        self.P_priori = np.zeros((9,9))
        self.F = np.zeros((9,9))
        self.H = np.zeros((4,9))
    
        self.current_time = 0
        self.s_dt = dt

        self.x_k = np.array([pos_x, vel_x, accel_x, pos_y, vel_y, accel_y, pos_z, vel_z, accel_z]).T

        for i in range(3):
            # self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
            #                                 [0.0, 1.0, dt],
            #                                 [0.0, 0.0, 1.0]]
            # Q = np.array([[(dt**5)/20., ((dt**4)/8.)*80., (dt**3)/6.],
            #               [((dt**4)/8.)*80., ((dt**3)/8.), (dt**2)/2.],
            #               [(dt**3)/6., (dt**2)/2., dt]])
            # self.Q[3*i:3*i+3, 3*i:3*i+3] = Q*13.
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)


        # accelerometer 3 axes
        # barometric altimeter 1 axis
        self.H = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

    def priori(self, R: np.ndarray, T:float, m:float, r:float, h:float, Cn:float, Ca:float, Cp:float, rho:float, bno_attitude:tuple, accel:tuple, timestep):
        """Sets priori state and covariance
            Try reading this: https://en.wikipedia.org/wiki/Kalman_filter#Details
            But basically predicts step
        Args:
            u (float): control input
        """
        # states tracked: x, vx, ax, y, vy, ay, z, vz, az
        # pos in x_k is every third element starting from 0

        pos_x, pos_y, pos_z = self.x_k[0], self.x_k[3], self.x_k[6]
        vel_x, vel_y, vel_z = self.x_k[1], self.x_k[4], self.x_k[7]
        vel_mag = np.linalg.norm([vel_x, vel_y, vel_z])
        acc_x, acc_y, acc_z = accel
        w_acc = vct.world_to_body(*bno_attitude, np.array([acc_x, acc_y, acc_z])) + np.array([-9.81, 0, 0])

        # approximate angular velocity
        w_x, w_y, w_z = w_acc * self.s_dt
        J_x = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2

        J_z = J_y

        Fax, Fay, Faz = 0,0,0                                   # aerodynamic forces expressed on the body in each direction
        Fax = -0.5*rho*(vel_mag**2)*float(Ca)*(np.pi*r**2)             # drag force

        #TODO: Verify Cn is being pulled from Aneesh's lookup table
        Fay = 0.5*rho*(vel_mag**2)*Cn*(np.pi*r**2) # take a look at only the Cn and see what happens
        Faz = Fay

        g = 9.81 # Earth gravity
        Fg = np.array([-g, 0, 0])
        Fg_body = np.linalg.inv(R) @ Fg
        Fgx, Fgy, Fgz = Fg_body[0], Fg_body[1], Fg_body[2]
        Ftx, Fty, Ftz = T[0],T[1],T[2]

        # states tracked: x, vx, ax, y, vy, ay, z, vz, az
        xdot = np.array([vel_x, (Fax + Ftx + Fgx) / m - (w_y*vel_z - w_z*vel_y), 0.0,
                 vel_y, (Fay + Fty + Fgy) / m - (w_z*vel_x - w_x*vel_z), 0.0,
                 vel_z, (Faz + Ftz + Fgz) / m - (w_x*vel_y - w_y*vel_x), 0.0
                ])
        self.x_priori = self.x_k + xdot * self.s_dt
        self.F = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0], 
                           [0, 0, 0, 0, w_z, 0, 0, -w_y, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                           [0, 0, 0, 0, 1, 0, 0, 0, 0], 
                           [0, -w_z, 0, 0, 0, 0, 0, w_x, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 1, 0], 
                           [0, w_y, 0, 0, -w_x, 0, 0, 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0]])

        self.P_priori = self.F @ self.P_k @ self.F.T + self.Q

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
        w_acc = vct.body_to_world(*bno_attitude, np.array([x_accel, y_accel, z_accel])) + np.array([-9.81, 0, 0])
        y_k = np.array([x_pos, *acc]).T

        self.x_k = self.x_priori + K @ (y_k - (self.H @ self.x_priori))
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
        self.x_k[1:3] = [0, 0]
        


