# ______      _____ _              __________ ___________ 
# | ___ \    /  ___(_)            / ___|  _  \  _  |  ___|
# | |_/ /   _\ `--. _ _ __ ___   / /___| | | | | | | |_   
# |  __/ | | |`--. \ | '_ ` _ \  | ___ \ | | | | | |  _|  
# | |  | |_| /\__/ / | | | | | | | \_/ | |/ /\ \_/ / |    
# \_|   \__, \____/|_|_| |_| |_| \_____/___/  \___/\_|    
#        __/ |                                            
#       |___/                                             

# A 6DOF RK-4 Based simulation that uses RASAero aerodynamic data and known motor thrust data to simulate motion of the rocket
# with simulated sensor data as well as a implementation of the Extended Kalman Filter and active drag PID controller for the ISS
# Spaceshot entry for the 2023 IREC competition. This simulation was used to quantify the effects of the airbrakes, test 
# different system design methodlogies, and provide preliminary tuning for the EKF and controller prior to implementation in 
# SILSIM and flight software.
# 
# 2022-2023 Guidance, Navigation, and Control Main Contributors #
# Sub-Team Lead: Parth Shrotri (2024)
# Colin Kinsey (2024)
# Evan Yu (2025)
# Rithvik Bhogavilli (2025)
# Kabir Cheema (2025)
# Freya Bansal (2025)
# Ishaan Bansal (2025)
# Ethan Pereira (2026)

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import shutil

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import estimation.ekf as ekf
import estimation.r_ekf as r_ekf
import properties.properties as prop
import properties.data_loader as dataloader
import simulator as sim_class
import dynamics.sensors as sensors
import time
import estimation.apogee_estimator as apg
import dynamics.rocket as rocket_model
import environment.atmosphere as atmosphere

# Load desired config file
config = dataloader.config

# Runs simulation for the specific component of the rocket

gl_apogee_last = 0
gl_apogee_d1_last = 0
gl_apogee_d2_last = []
gl_apogee_plot = []
gl_apogee_plot_d1 = []
gl_apogee_plot_d2 = []
gl_apogee_d2_accum = []
gl_apogee_accum_plot = []

class Simulation:
    # dt can be dynamic in the future, so we need to 
    def __init__(self, rocket: rocket_model.Rocket, motor, dt, x0, time_stamp=0, stages=[]):
        self.rocket = rocket
        self.stages = stages
        self.motor = motor
        self.dt = dt
        self.x = x0.copy()
        self.baro = 0
        self.time_stamp = time_stamp
        self.sensor_config = self.rocket.stage_config['sensors']
        self.init_kalman_filters()

    def init_kalman_filters(self):
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
        self.apogee_estimator = apg.Apogee(self.kalman_filter.get_state(), 0.1, 0.01, 3, 30, atm, self.rocket)

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

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, 0, current_state[0], self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, dt)
            self.time_step()

    def execute_stage(self):
        global gl_apogee_accum_plot, gl_apogee_plot, gl_apogee_plot_d2, gl_apogee_plot_d1
        # Run the stages
        stage_separation_delay = 1

        minimum_ignition_delay = 10
        motor_ignition_delay = 0
        maximum_ignition_delay = 15
        is_ignited = False

        stage_time = self.time_stamp

        time_stamps = []

        start = True
        print(f"Staged at {self.time_stamp}")
        apogee_optimizer = apg.ApogeeOptimizer(1,1,1,1)
        while self.time_stamp < stage_time + self.rocket.get_motor().get_burn_time() + stage_separation_delay + motor_ignition_delay:
            if start:
                start = False


            # Get sensor data
            baro_alt, accel, gyro, bno_ang_pos = self.get_sensor_data()
            self.update_kalman(baro_alt, accel, gyro, bno_ang_pos)
            current_state, current_covariance, current_state_r = self.get_kalman_state()

            # print(self.rocket.current_stage)
            # Apogee optimization
            if self.rocket.current_stage == 0 and not is_ignited:
                
                # Only do this for sustainer stage

                # Get apogee estimate:
                apogee_estimate = self.apogee_estimator.predict_apogee(current_state[0:3])
                apogee_optimizer.add_sample(apogee_estimate, dt)

                print(f"Coasting: {(self.time_stamp - stage_time):.2f}s | Current motor delay: {motor_ignition_delay:.2f}s | Current Apogee Estimate: {apogee_estimate:.2f}m", end="\r")

                gl_apogee_accum_plot.append(apogee_optimizer.get_d2_smooth())
                gl_apogee_plot.append(apogee_estimate)
                gl_apogee_plot_d2.append(apogee_optimizer.get_d2())
                gl_apogee_plot_d1.append(apogee_optimizer.get_d1())
                time_stamps.append(self.time_stamp)

                if motor_ignition_delay < minimum_ignition_delay:
                    motor_ignition_delay = minimum_ignition_delay

                if motor_ignition_delay > maximum_ignition_delay:
                    motor_ignition_delay = maximum_ignition_delay

                if self.time_stamp > motor_ignition_delay + stage_time:

                    plt.plot(time_stamps, gl_apogee_plot, label="apogee_estimate")
                    plt.plot(time_stamps, gl_apogee_plot_d1, label="apogee_d1")
                    plt.plot(time_stamps, gl_apogee_plot_d2, label="apogee_d2")
                    plt.plot(time_stamps, gl_apogee_accum_plot, label="apogee_smooth_d2")
                    plt.legend()

                    plt.show()
                    print("")
                    print(f"(delayed) Motor ignited at {self.time_stamp}")
                    is_ignited = True
                    self.rocket.get_motor().ignite(self.time_stamp)
                    ignition_time = self.time_stamp
            else:
                if not is_ignited:
                    # If not second stage, immediately ignite.
                    print(f"Motor ignited at {self.time_stamp}")
                    is_ignited = True
                    self.rocket.get_motor().ignite(self.time_stamp)
                    ignition_time = self.time_stamp



            apogee_est = 0

            self.rocket.set_motor_mass(self.time_stamp)

            self.x, alpha = sim.RK4(self.x, dt, self.time_stamp, 0)

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, alpha, apogee_est, self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, dt)
            self.time_step()


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

            # apogee_est = self.apogee_estimator.predict_apogee(current_state[0:3])

            self.x, alpha = sim.RK4(self.x, dt, self.time_stamp, 0)

            self.rocket.add_to_dict(self.x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_covariance, current_state_r, alpha, apogee_est, self.rocket.get_rocket_dry_mass(), self.rocket.get_total_motor_mass(self.time_stamp), 0, dt)
            self.time_step()
    
def simulator(x0, rocket, motor, dt) -> None:
    '''Method which handles running the simulation and logging sim data to dict

    Args:
        x0 (np.array): state vector initialized to 0s [6x3]
             x:         y:         z:
           [[pos,       pos,       pos],
            [vel,       vel,       vel],
            [accel,     accel,     accel],
            [ang_pos,   ang_pos,   ang_pos],
            [ang_vel,   ang_vel,   ang_vel],
            [ang_accel, ang_accel, ang_accel]]
        dt (float): time step between each iteration in simulation
    '''
    simulator = Simulation(rocket, motor, dt, x0, stages=rocket.stages)
    simulator.idle_stage()
    t_start = time.time()
    simulator.run_stages()
    simulator.coast()
    t_end = time.time() - t_start
    print(f"Runtime: {t_end:.2f} seconds")


if __name__ == '__main__':
    x0 = np.zeros((6, 3))
    x0[3] = [0, 0.05, 0]
    dt = 0.01

    atm = atmosphere.Atmosphere(enable_direction_variance=True, enable_magnitude_variance=True)

    stages = []
    for stage in config['rocket']['stages'][1:]:
        stages.append(rocket_model.Rocket(stage, atm=atm))
    rocket = rocket_model.Rocket(config['rocket']['stages'][0], atm=atm, stages=stages)

    motor = rocket.motor
    sim = sim_class.Simulator(atm=atm, rocket=rocket)

    simulator(x0, rocket, motor, dt)

    print("Writing to file...")

    record = rocket.to_csv()
    output_file = os.path.join(os.path.dirname(__file__), config["meta"]["output_file"])
    with open(output_file, 'w') as f:
        f.write("time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,accel_x,accel_y,accel_z,ang_pos_x,ang_pos_y,ang_pos_z,ang_vel_x,ang_vel_y,ang_vel_z,ang_accel_x,ang_accel_y,ang_accel_z,alpha,rocket_total_mass,motor_mass,flap_ext,baro_alt,imu_accel_x,imu_accel_y,imu_accel_z,imu_ang_pos_x,imu_ang_pos_y,imu_ang_pos_z,imu_gyro_x,imu_gyro_y,imu_gyro_z,apogee_estimate,kalman_pos_x,kalman_vel_x,kalman_accel_x,kalman_pos_y,kalman_vel_y,kalman_accel_y,kalman_pos_z,kalman_vel_z,kalman_accel_z,pos_cov_x,vel_cov_x,accel_cov_x,pos_cov_y,vel_cov_y,accel_cov_y,pos_cov_z,vel_cov_z,accel_cov_z,kalman_rpos_x,kalman_rvel_x,kalman_raccel_x,kalman_rpos_y,kalman_rvel_y,kalman_raccel_y,kalman_rpos_z,kalman_rvel_z,kalman_raccel_z\n")
        for point in record:
            f.write(f"{','.join(point)}\n")
