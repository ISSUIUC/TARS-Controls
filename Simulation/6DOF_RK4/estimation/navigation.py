# import ekf, rekf
import sys
import os
import numpy as np
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
import estimation.ekf as ekf
import estimation.r_ekf as r_ekf
import dynamics.sensors as sensors

# create instances of kalman filter classes in ekf and rekf
# take the rocket state tracking from main/simulation to rocket, rocket itself to have a state instead of the simulation
# rocket has instance of the navigation class
class Navigation:
    def __init__(self, dt, sensor_config,x0):
        # TODO: Init with previous rocket data
        self.kalman_filter = ekf.KalmanFilter(dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.sensor_config = sensor_config
        #self.ekf = ekf
        #self.r_ekf = r_ekf
        self.x = x0
        # TODO: init this data with the previous rocket becuase of staging
        x_data = []
        y_data = []
        z_data = []
        for i in range(10):
            reading = sensors.get_accelerometer_data(self.x, self.sensor_config)
            x_data.append(reading[0])
            y_data.append(reading[1])
            z_data.append(reading[2])

        # accel_tracker = np.array([])

        ax = sum(x_data)/len(x_data)
        ay = sum(y_data)/len(y_data)
        az = sum(z_data)/len(z_data)

        pitch = -1 * (np.arctan2(-az,-ay) + np.pi/2)
        yaw = np.arctan2(-ax,-ay) + np.pi/2
        self.r_kalman_filter = r_ekf.KalmanFilter_R(dt, 0.0, 0.0, 0.0, pitch, 0.0, 0.0, yaw, 0.0, 0.0)
        ## self.apogee_estimator = apg.Apogee(self.kalman_filter.get_state(), 0.1, 0.01, 3, 30, atm, self.rocket.stage_config)
    def update_state(self,baro_alt, accel, gyro, bno_ang_pos):
        # x_pos, vel_x,  vel_y,  vel_z, alpha_x, alpha_y, alpha_z
        bno_attitude = bno_ang_pos
        x_pos = baro_alt    
        x_accel, y_accel, z_accel = accel[0], accel[1], accel[2]
        self.kalman_filter.priori()
        self.r_kalman_filter.priori()
        self.kalman_filter.update(bno_attitude, x_pos, x_accel, y_accel, z_accel)
        self.r_kalman_filter.update(*gyro, *accel)
