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
        
        """ for writing change in rotation and shi
        Aroll, Apitch, Ayaw : Aerodynamic moment vector for roll, pitch, yaw
        Proll, Ppitch, Pyaw : Propulsion moment vector for roll, pitch, yaw 
        Ix , Iy, Iz         : Compon. of inertia (diagonal elements of inertia matrix when inertia products are zero)
        0roll, 0pitch, 0yaw : angular rate of change for roll, pitch, yaw


        roll_ddot = Aroll + Proll - (pitch*yaw)*(Iz - Iy)
                    -------------------------------------
                                    Ix
        
        pitch_ddot = Apitch + Ppitch - (roll*yaw)*(Ix - Iz)
                    -------------------------------------
                                    Iy

        yaw_ddot = Ayaw + Pyaw - (roll*pitch)*(Iy - Ix)
                    -------------------------------------
                                    Iz

        Aroll = 0.5*rho*V^2*Cl*A
        Apitch = 0.5*rho*V^2*Cm*A
        Ayaw = 0.5*rho*V^2*Cn*A

        Cl = roll moment coefficient
        Cm = pitch moment coefficient
        Cn = yaw moment coefficient

        Cl = C_ls S_r + d/2V_m (C_lp P)
        Cm = C_mref - C_Nz (x_cm - x_ref)/d + d/2V_m(C_mq + C_malpha)q dimensionless
        Cn = C_nref + C_Ny (x_cm - x_ref)/d + d/2V_m(C_nr + C_nbeta)r

        C_lp = roll damping derivative relative to roll rate p, rad^-1 (deg^-1)
        C_ls = slope of curve formed by roll moment coefficient C1 versus control surface deflection, rad^-1 (deg^-1)
        C_mref = pitching moment coefficient about reference moment station, dimensionless.
        C_mq = pitch damping derivatives relative to pitch rate q, rad^-1 (deg^-1)
        C_malpha = pitch damping derivative relative to angle of attack rate α̇ 
        (slope of curve formed by C_alpha verses alpha), rad^-1 (deg^-1)
        C_Ny = coefficient corresponding to component of normal force on yb-axis, dimensionless.
        C_Nz = coefficient corresponding to component of normal force on zb-axis, dimensionless.
        C_n = yaw damping derivative relative to yaw rate r_dot, rad^-1 (deg^-1)
        C_nref = yawing moment coefficient about reference moment station, dimensionless.
        C_nbeta = yaw damping derivative relative to angle of sideslip rate beta_dot, rad^-1 (deg^-1)
        x_cm = instantaneous distance from rocket nose to center of mass, m
        x_ref = distance from rocket nose to reference moment station, m
        S_r = effective control surface deflection causing rolling moment, rad(deg)

        """
        


    # TODO: Fix priori and update step, make sure that names make sense, and add comments.
    def priori(self, m:float, r:float, h:float, rho:float, cm:float, cp:float, roll_p:float, pitch_p:float, yaw_p:float, x_rot):
        J_x = 1/2 * m * r**2
        J_y = 1/3 * m * h**2 + 1/4 * m * r**2

        J_z = J_y
        
        
        Croll_a = 1 #aerodynamic roll moment coefficient
        Cpitch_a = 1 #aerodynamic pitch moment coefficient
        Cyaw_a = 1 #aerodynamic yaw moment coefficient

        vel_x, vel_y, vel_z = self.x_k[1], self.x_k[4], self.x_k[7]
        vel_mag = np.linalg.norm([vel_x, vel_y, vel_z])

        roll_a = 0.5*rho*vel_mag*Croll_a*np.pi*(r)**2
        pitch_a = 0.5*rho*vel_mag*Cpitch_a*np.pi*(r)**2
        yaw_a = 0.5*rho*vel_mag*Cyaw_a*np.pi*(r)**2

        pos_roll = self.x_k[0]
        pos_pitch = self.x_k[1]
        pos_yaw = self.x_k[2]

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

            # vel_mag = (vel_x**2 + vel_y**2 + vel_z**2)**0.5

            # For thrust calculations, mutiply Thrust by a rotation matrix in rollpitchyaw directions
            # then take the cross product of moment arm with thrust to find the forces in relative directions
            # To find moment arm distance, take scalar diff of cm and cp, and make a vector [x, 0, 0], then multiply by rotation matrix
            # finally take the cross product between the two and it will return roll pitch yaw thrust moment coefficients
        ])




        F = np.array([[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1],
                        [0, 0, 0, 0, -w_z*(-J_y + J_z)/J_x, -w_y*(-J_y + J_z)/J_x], [0, 0, 0, -w_z*(J_x - J_z)/J_y, 0, -w_x*(J_x - J_z)/J_y], [0, 0, 0, -w_y*(-J_x + J_y)/J_z, -w_x*(-J_x + J_y)/J_z, 0]])
        self.x_priori = self.F @ self.x_k
        self.P_priori = (self.F @ self.P_k @ self.F.T) + self.Q

    def skew(self, v):
        return np.array([[0, -v[2], v[1]],
                        [v[2], 0, -v[0]],
                        [-v[1], v[0], 0]])

    def update(self, vel_x,  vel_y,  vel_z, x_accel, y_accel, z_accel):
        K = (self.P_priori @ self.H.T) @ np.linalg.inv(self.H @ self.P_priori @ self.H.T + self.R)
        body_rot_rate = np.array([vel_x,vel_y,vel_z])
        world_rot_rate = vct.body_to_world(self.x_k[0], self.x_k[3], self.x_k[6], body_rot_rate)
        roll = 0
        yaw = np.arctan2(y_accel,x_accel)
        pitch = np.arctan2(z_accel, np.sqrt(y_accel**2 + x_accel**2))
        #Euler-Poisson
        o_b_to_a = vct.body_to_world(roll, pitch, yaw, [1,0,0])
        o_dot = (-self.skew(body_rot_rate)).T @ o_b_to_a
        current_euler_orientation = (o_dot * self.dt) + o_b_to_a
        y_k = np.array([0, world_rot_rate[0], pitch, world_rot_rate[1], yaw, world_rot_rate[2]]).T

        

        self.x_k = self.x_priori + K @ (y_k - self.H @ self.x_priori)
        self.P_k = (np.eye(len(K)) - K @ self.H) @ self.P_priori

        self.current_time += self.s_dt

    def get_state(self):
        return self.x_k
