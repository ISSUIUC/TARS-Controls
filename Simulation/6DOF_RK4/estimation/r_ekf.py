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
        self.P_k = np.eye(9)
        self.x_priori = np.zeros((9,1))
        self.P_priori = np.zeros((9,9))
        self.F = np.zeros((9,9))
        self.H = np.zeros((6,9))

        # used for aerodynamic force calc
        self.vel = np.zeros(3)

        self.current_time = 0
        self.s_dt = dt


        
        self.x_k = np.array([[roll, w_x, a_x, 
                             pitch, w_y, a_y, 
                             yaw, w_z, a_z]]).T

        for i in range(3):
            # why is F initialized here? didnt see this in the onb-errorekfraa branch for ekf.py
            self.F[3*i:3*i+3, 3*i:3*i+3] = [[1.0, dt, (dt**2) / 2],
                                            [0.0, 1.0, dt],
                                            [0.0, 0.0, 1.0]]
            self.Q[3*i:3*i+3, 3*i:3*i+3] = Q_continuous_white_noise(3, dt, 13.)
            

        # double check if this is right
        self.H = np.array([
            [1, 0, 0, 0, 0, 0, 0, 0, 0],  
            [0, 1, 0, 0, 0, 0, 0, 0, 0],  
            [0, 0, 0, 1, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 1, 0, 0, 0, 0],  
            [0, 0, 0, 0, 0, 0, 1, 0, 0],  
            [0, 0, 0, 0, 0, 0, 0, 1, 0],  
        ])


        
    
        


    # TODO: Fix priori and update step, make sure that names make sense, and add comments.
    def priori(self, R, m, T, r, h, rho, cm, cp, Cx_aero, Cy_aero, Cz_aero, bno_attitude, x_accel, y_accel, z_accel):
        J_x = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2
        J_z = J_y

        grav = -9.81
        g = np.array([grav, 0, 0])
        gBody = vct.world_to_body(*bno_attitude, g)

        cp_val = float(cp)
        thrust = T    
        # this jawn was printing out to be 128.45, 0, 0. guessing its cm so i divided by 100
        cP = np.array([cp_val/100, 0.0 , 0.0])
        momentVector =  cP - cm

        # acc = vct.body_to_world(*bno_attitude, np.array([x_accel, y_accel, z_accel])) + np.array([-9.81, 0, 0])
        accBody = np.array([x_accel,y_accel,z_accel])
        aBody = accBody - gBody
        self.vel += aBody * self.dt
        vel_mag = np.linalg.norm(self.vel)

        self.x_k = self.x_k.flatten()
        w_x, w_y, w_z = self.x_k[1].item(), self.x_k[4].item(), self.x_k[7].item()


        A = np.pi * r**2
        # x_force = -0.5 * rho * self.vel[0]**2 * Cx_aero * A * np.sign(self.vel[0])
        # y_force = -0.5 * rho * self.vel[1]**2 * Cy_aero * A * np.sign(self.vel[1])
        # z_force = -0.5 * rho * self.vel[2]**2 * Cz_aero * A * np.sign(self.vel[2])

        x_force = -0.5 * rho * vel_mag * Cx_aero * A * np.sign(self.vel[0])
        y_force = -0.5 * rho * vel_mag * Cy_aero * A * np.sign(self.vel[1])
        z_force = -0.5 * rho * vel_mag * Cz_aero * A * np.sign(self.vel[2])


        aForce = np.array([x_force, y_force, z_force])  # body frame
        aMoment = np.cross(momentVector, aForce)
        x_aero, y_aero, z_aero = aMoment

        Mt = np.cross(momentVector, thrust)  
        Mtx = Mt[0]; Mty = Mt[1]; Mtz = Mt[2]
            
        # atleast this is setup right 💀
        xdot = np.array([
            w_x, 
            (x_aero + Mtx - w_y*w_z*(J_z - J_y)) / J_x, 
            0.0, 
            w_y, 
            (y_aero + Mty - w_x*w_z*(J_x - J_z)) / J_y, 
            0.0,
            w_z, 
            (z_aero + Mtz - w_x*w_y*(J_y - J_x)) / J_z, 
            0.0 
        ])

        

        self.F = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, -w_z*(-J_y + J_z)/J_x, 0, 0, -w_y*(-J_y + J_z)/J_x, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 1, 0, 0, 0, 0], 
                      [0, -w_z*(J_x - J_z)/J_y, 0, 0, 0, 0, 0, -w_x*(J_x - J_z)/J_y, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 1, 0], 
                      [0, -w_y*(-J_x + J_y)/J_z, 0, 0, -w_x*(-J_x + J_y)/J_z, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0, 0, 0]])

        #self.x_priori = self.F @ self.x_k
        self.x_priori = self.x_k + xdot * self.s_dt
        self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q

    def update(self, vel_x,  vel_y,  vel_z, x_accel, y_accel, z_accel, bno_roll, bno_pitch, bno_yaw):
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)

        body_rot_rate = np.array([vel_x,vel_y,vel_z])
        world_rot_rate = vct.body_to_world(self.x_k[0], self.x_k[3], self.x_k[6], body_rot_rate)
        # yaw = np.arctan2(y_accel,x_accel)
        # pitch = np.arctan2(z_accel, np.sqrt(y_accel**2 + x_accel**2))
        # y_k = np.array([bno_roll, world_rot_rate[0], bno_pitch, world_rot_rate[1], bno_yaw, world_rot_rate[2]]).T
        y_k = np.array([bno_roll, vel_x, bno_pitch, vel_y, 0, vel_z]).T

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori 
        self.current_time += self.s_dt


    def get_state(self):
        return self.x_k


