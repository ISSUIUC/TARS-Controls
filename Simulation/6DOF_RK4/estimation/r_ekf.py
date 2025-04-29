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


        
        self.x_k = np.array([[roll, w_x, a_x, 
                             pitch, w_y, a_y, 
                             yaw, w_z, a_z]]).T

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
    def priori(self, R: np.ndarray, m:float, T:float, r:float, h:float, rho:float, cm:float, cp:float, Cx_aero: float, Cy_aero: float, Cz_aero: float):
        # m, r, and h are being returned as floats
        J_x = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2
        J_z = J_y

        # mass, height, radius are all floats

        # coming up as floats so we good
        w_x, w_y, w_z = self.x_k[1].item(), self.x_k[4].item(), self.x_k[7].item()
        w_mag = np.linalg.norm([w_x, w_y, w_z])

        if Cy_aero > 0:
            print("cyaero: ", Cy_aero)

        x_aero = 0.5 * rho * w_mag * Cx_aero * np.pi*(r)**2
        y_aero = 0.5 * rho * w_mag * Cy_aero * np.pi*(r)**2
        z_aero = 0.5 * rho * w_mag * Cz_aero * np.pi*(r)**2

        cp_val = float(cp)
        thrust = T    
        cP = np.array([cp_val, 0.0 , 0.0])
        vec_cp_to_cm =  cP - cm

        thrust_world = R @ thrust
        vec_cp_cm_world = R @ vec_cp_to_cm

        Mt = np.cross(vec_cp_cm_world, thrust_world)  
        Mtx = Mt[0]; Mty = Mt[1]; Mtz = Mt[2]

        xdot = np.array([
            [w_x], 
            [(x_aero + Mtx - w_y*w_z*(J_z - J_y)) / J_x], 
            [1.0], 
            [w_y], 
            [(y_aero + Mty - w_x*w_z*(J_x - J_z)) / J_y], 
            [1.0],
            [w_z], 
            [(z_aero + Mtz - w_x*w_y*(J_y - J_x)) / J_z], 
            [1.0] 
        ])
        if xdot[1].item() > 1:
            print("roll position: ", xdot[1])

        self.F = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, -w_z*(-J_y + J_z)/J_x, 0, 0, -w_y*(-J_y + J_z)/J_x, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 1, 0, 0, 0, 0], 
                      [0, -w_z*(J_x - J_z)/J_y, 0, 0, 0, 0, 0, -w_x*(J_x - J_z)/J_y, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 1, 0], 
                      [0, -w_y*(-J_x + J_y)/J_z, 0, 0, -w_x*(-J_x + J_y)/J_z, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0]])

        self.x_priori = self.x_k + xdot * self.s_dt
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


