import sys
import os
import numpy as np
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
import estimation.ekf as ekf
import estimation.r_ekf as r_ekf
import dynamics.sensors as sensors
## Class for the Rocket class and object to interact with the Kalman filter (and eventually staging optimization)
## Pysim and other files should not directly interact with Navigation or the ekf/r_ekf files but should instead utilize this class through a rocket object
class Navigation:
    def __init__(self, dt, sensor_config,x0):
        # TODO: Init with previous rocket data
        self.kalman_filter = ekf.KalmanFilter(dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.sensor_config = sensor_config
        self.x = x0
        # TODO: init this data with the previous rocket becuase of staging 
        x_data = np.empty(10)
        y_data = np.empty(10)
        z_data = np.empty(10)
        for i in range(10):
            reading = sensors.get_accelerometer_data(self.x, self.sensor_config)
            x_data[i] = reading[0]
            y_data[i] = reading[1]
            z_data[i] = reading[2]
        ## Acceleration value given by an average of a sum of values
        ax = sum(x_data)/len(x_data)
        ay = sum(y_data)/len(y_data)
        az = sum(z_data)/len(z_data)
        ## Angular component calculations
        pitch = -1 * (np.arctan2(-az,-ay) + np.pi/2)
        yaw = np.arctan2(-ax,-ay) + np.pi/2
        self.r_kalman_filter = r_ekf.KalmanFilter_R(dt, 0.0, 0.0, 0.0, pitch, 0.0, 0.0, yaw, 0.0, 0.0)
        ## self.apogee_estimator = apg.Apogee(self.kalman_filter.get_state(), 0.1, 0.01, 3, 30, atm, self.rocket.stage_config)
        
    def update_state(self,baro_alt, accel, gyro, bno_ang_pos):
        self.kalman_filter.priori()
        self.r_kalman_filter.priori()
        self.kalman_filter.update(bno_ang_pos, baro_alt, accel[0], accel[1], accel[2])
        self.r_kalman_filter.update(*gyro, *accel)
