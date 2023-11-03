import numpy as np

import estimation.ekf as ekf
import estimation.r_ekf as r_ekf
import dynamics.sensors as sensors
import estimation.apogee_estimator as apg

class Simulation:
    # dt can be dynamic in the future, so we need to 
    def __init__(self, sim, rocket, atm, dt, x0, kf_dt, time_stamp=0, stages=[]):
        self.sim = sim
        self.rocket = rocket
        self.stages = stages
        self.motor = rocket.motor
        self.atm = atm
        self.dt = dt
        self.x = x0.copy()
        self.baro = 0
        self.time_stamp = time_stamp
        self.sensor_config = self.rocket.stage_config['sensors']
        self.init_gnc(kf_dt)

    def init_gnc(self, dt):
        # TODO: Init with previous rocket data
        self.kalman_filter = ekf.KalmanFilter(dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        # TODO: init this data with the previous rocket becuase of staging
        x_data = []
        y_data = []
        z_data = []
        for i in range(10):
            reading = sensors.get_accelerometer_data(self.x, self.sensor_config)
            x_data.append(reading[0])
            y_data.append(reading[1])
            z_data.append(reading[2])

        accel_tracker = np.array([])

        ax = sum(x_data)/len(x_data)
        ay = sum(y_data)/len(y_data)
        az = sum(z_data)/len(z_data)

        pitch = -1 * (np.arctan2(-az,-ay) + np.pi/2)
        yaw = np.arctan2(-ax,-ay) + np.pi/2
        self.r_kalman_filter = r_ekf.KalmanFilter_R(dt, 0.0, 0.0, 0.0, pitch, 0.0, 0.0, yaw, 0.0, 0.0)
        self.apogee_estimator = apg.Apogee(self.kalman_filter.get_state(), 0.1, 0.01, 3, 30, self.atm, self.rocket.stage_config)

    def update_kalman(self, baro_alt, accel, gyro, bno_ang_pos):
        self.kalman_filter.priori()
        self.kalman_filter.update(bno_ang_pos, baro_alt,
                            accel[0], accel[1], accel[2])
        
        self.r_kalman_filter.priori()
        self.r_kalman_filter.update(*gyro, *accel)
    
    def time_step(self):
        self.time_stamp += self.dt

    def idle_stage(self):
        while self.time_stamp < self.rocket.delay:
            baro_alt, accel, gyro, bno_ang_pos = self.get_sensor_data()
            self.update_kalman(baro_alt, accel, gyro, bno_ang_pos)
            self.kalman_filter.reset_lateral_pos()
            current_state, current_covariance, current_state_r = self.get_kalman_state()

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, 0, current_state[0], self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, self.dt)
            self.time_step()

    def execute_stage(self, print=False):
        # Run the stages
        stage_separation_delay = 1
        self.rocket.get_motor().ignite(self.time_stamp)

        ignition_time = self.time_stamp
        start = True
        if print:
            print(f"Staged at {self.time_stamp}")
        while self.time_stamp < ignition_time + self.rocket.get_motor().get_burn_time() + stage_separation_delay:
            # Get sensor data
            baro_alt, accel, gyro, bno_ang_pos = self.get_sensor_data()
            self.update_kalman(baro_alt, accel, gyro, bno_ang_pos)
            current_state, current_covariance, current_state_r = self.get_kalman_state()

            apogee_est = self.apogee_estimator.predict_apogee(current_state[0:3])

            self.rocket.set_motor_mass(self.time_stamp)

            is_staging = start and self.rocket.current_stage != -1
            self.x, alpha = self.sim.RK4(self.x, self.dt, self.time_stamp, is_staging, 0)

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, alpha, apogee_est, self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, self.dt)
            self.time_step()
            if start:
                start = False

    def run_stages(self):
        has_more_stages = True
        while has_more_stages:
            self.execute_stage()
            has_more_stages = self.rocket.separate_stage(self.time_stamp)

    # Function to retrive all sensor data
    def get_sensor_data(self):
        return (sensors.get_barometer_data(self.x, self.sensor_config),
                sensors.get_accelerometer_data(self.x, self.sensor_config),
                sensors.get_gyro_data(self.x, self.sensor_config), 
                sensors.get_bno_orientation(self.x, self.sensor_config))

    def get_kalman_state(self):
        current_state = self.kalman_filter.get_state()
        current_cov = self.kalman_filter.get_covariance()
        current_state_r = self.r_kalman_filter.get_state()
        return (current_state, current_cov, current_state_r)

    def coast(self):
        while self.x[1, 0] >= 0:
        # Get sensor data
            baro_alt, accel, gyro, bno_ang_pos = self.get_sensor_data()
            self.update_kalman(baro_alt, accel, gyro, bno_ang_pos)
            current_state, current_covariance, current_state_r = self.get_kalman_state()

            apogee_est = self.apogee_estimator.predict_apogee(current_state[0:3])

            self.x, alpha = self.sim.RK4(self.x, self.dt, self.time_stamp, 0)

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, alpha, apogee_est, self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, self.dt)
            self.time_step()
    