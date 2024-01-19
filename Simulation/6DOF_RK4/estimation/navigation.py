# import ekf, rekf
import sys
import os
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
import ekf
import r_ekf
# create instances of kalman filter classes in ekf and rekf
# take the rocket state tracking from main/simulation to rocket, rocket itself to have a state instead of the simulation
# rocket has instance of the navigation class
class Navigation:
    def __init__(self, dt, p_x, v_x, a_x, p_y, v_y, a_y, p_z, v_z, a_z, theta_x, omega_x, alpha_x, theta_y, omega_y, alpha_y, theta_z, omega_z, alpha_z):
        self.ekf = ekf.KalmanFilter(dt, p_x, v_x, a_x, p_y, v_y, a_y, p_z, v_z, a_z)
        self.r_ekf = r_ekf.KalmanFilter(dt, theta_x, omega_x, alpha_x, theta_y, omega_y, alpha_y, theta_z, omega_z, alpha_z)
    def update_state(self,bno_attitude, x_pos, x_accel, y_accel, z_accel, vel_x,  vel_y,  vel_z, alpha_x, alpha_y, alpha_z):
        self.ekf.priori()
        self.r_ekf.priori()
        self.ekf.update(bno_attitude, x_pos, x_accel, y_accel, z_accel)
        self.r_ekf.update(vel_x,  vel_y,  vel_z, alpha_x, alpha_y, alpha_z)