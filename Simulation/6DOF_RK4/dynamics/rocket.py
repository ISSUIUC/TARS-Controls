import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import dynamics.motor as motor
import dynamics.forces as forces

class Rocket:
    motor = None
    forces = None
    stage_config = None

    def init_dicts(self):
        self.sim_dict = {
            "pos": [],
            "vel": [],
            "accel": [],
            "ang_pos": [],
            "ang_vel": [],
            "ang_accel": [],
            "alpha": [],
            "flap_ext": [],
            "rocket_total_mass": [],
            "motor_mass": [],
            "time": [],
        }

        self.kalman_dict = {
            "x": [],
            "y": [],
            "z": [],
            "cov_x": [],
            "cov_y": [],
            "cov_z": [],
            "rx": [],
            "ry": [],
            "rz": [],
            "time": []
        }

        self.sensor_dict = {
            "baro_alt": [],
            "imu_accel_x": [],
            "imu_accel_y": [],
            "imu_accel_z": [],
            "imu_ang_pos_x": [],
            "imu_ang_pos_y": [],
            "imu_ang_pos_z": [],
            "imu_gyro_x": [],
            "imu_gyro_y": [],
            "imu_gyro_z": [],
            "apogee_estimate": []
        }

    def __init__(self, stage_config, atm=None):
        self.stage_config = stage_config
        self.cm_rocket = stage_config["rocket_body"]["structure_cm"]
        self.cm_motor = stage_config["motor"]["cm"]
        self.cm = stage_config["rocket_body"]["combined_cm"]
        self.cp = stage_config["rocket_body"]["combined_cp"]

        self.impulse = stage_config["motor"]["impulse"]
        self.motor_mass = stage_config["motor"]["motor_mass"]
        self.delay = stage_config["motor"]["delay"]
        self.motor_lookup_file = stage_config["motor"]["motor_lookup_file"]

        self.rocket_dry_mass = stage_config["rocket_body"]["dry_mass"]
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        self.r_r = stage_config["rocket_body"]["radius"]
        self.l = stage_config["rocket_body"]["length"]
        self.A = math.pi * self.r_r ** 2
        self.A_s = 2 * self.r_r * self.l
        self.max_ext_length = stage_config["flaps"]["max_ext_length"]
        self.atm = atm

        self.init_dicts()
        # need to change motor and forces constructors to refer to properties from this class rather than properties file
        self.motor = motor.Motor(self.rocket_dry_mass, 
                                 self.cm, 
                                 self.cm_rocket, 
                                 self.cm_motor, 
                                 self.rocket_dry_mass, 
                                 impulse=self.impulse, 
                                 mass=self.motor_mass, 
                                 delay=self.delay, 
                                 lookup_file=self.motor_lookup_file)
        self.forces = forces.Forces(self.max_ext_length,
                                    self.cm,
                                    self.cp,
                                    self.A,
                                    self.A_s,
                                    self.rocket_dry_mass,
                                    self.motor,
                                    stage_config["rocket_body"]["rasaero_lookup_file"],
                                    self.atm)

    def set_motor_mass(self, timestamp):
        """Sets the mass of the motor at a given time
        
        Args:
            timestamp (float): Time in seconds
        """
        self.motor_mass = self.motor.get_mass(timestamp)
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        
    def I(self, total_mass): 
        """Returns the inertia matrix of the rocket
        
        Args:
            total_mass (float, optional): Total mass of the rocket. Defaults to rocket_total_mass.
        
        Returns:
            np.array: Inertia matrix of the rocket
        """
        return np.diag([(1/2) * total_mass * self.r_r**2,
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2),
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2)])

    
    def I_inv(self, total_mass): 
        """Returns the inverse of the inertia matrix of the rocket
        
        Args:
            total_mass (float, optional): Total mass of the rocket. Defaults to rocket_total_mass.
        """
        return np.diag([1/((1/2) * total_mass * self.r_r**2),
                                       1/((total_mass/12) * (self.l**2 + 3*self.r_r**2)),
                                       1/((total_mass/12) * (self.l**2 + 3*self.r_r**2))])
    
    # Appends data to rocket's storage
    def add_to_dict(self, x, baro_alt, accel, bno_ang_pos, gyro, kalman_filter, kf_cov, kalman_filter_r, alpha, apogee_estimation, rocket_total_mass, motor_mass, flap_ext, dt):
        # Append to sensor_dict
        self.sensor_dict["baro_alt"].append(baro_alt)
        self.sensor_dict["imu_accel_x"].append(accel[0])
        self.sensor_dict["imu_accel_y"].append(accel[1])
        self.sensor_dict["imu_accel_z"].append(accel[2])
        self.sensor_dict["imu_ang_pos_x"].append(bno_ang_pos[0])
        self.sensor_dict["imu_ang_pos_y"].append(bno_ang_pos[1])
        self.sensor_dict["imu_ang_pos_z"].append(bno_ang_pos[2])
        self.sensor_dict["imu_gyro_x"].append(gyro[0])
        self.sensor_dict["imu_gyro_y"].append(gyro[1])
        self.sensor_dict["imu_gyro_z"].append(gyro[2])
        self.sensor_dict["apogee_estimate"].append(apogee_estimation)

        self.kalman_dict["x"].append(kalman_filter[0:3])
        self.kalman_dict["y"].append(kalman_filter[3:6])
        self.kalman_dict["z"].append(kalman_filter[6:9])
        self.kalman_dict["cov_x"].append(kf_cov[0:3])
        self.kalman_dict["cov_y"].append(kf_cov[3:6])
        self.kalman_dict["cov_z"].append(kf_cov[6:9])
        self.kalman_dict["rx"].append(kalman_filter_r[0:3])
        self.kalman_dict["ry"].append(kalman_filter_r[3:6])
        self.kalman_dict["rz"].append(kalman_filter_r[6:9])

        # Update Simulator Log
        self.sim_dict["pos"].append(x[0])
        self.sim_dict["vel"].append(x[1])
        self.sim_dict["accel"].append(x[2])
        self.sim_dict["ang_pos"].append(x[3])
        self.sim_dict["ang_vel"].append(x[4])
        self.sim_dict["ang_accel"].append(x[5])
        self.sim_dict["time"].append(self.sim_dict["time"][-1] + dt if len(self.sim_dict["time"]) > 0 else 0)
        self.sim_dict["flap_ext"].append(flap_ext)                    
        self.sim_dict["alpha"].append(alpha)
        self.sim_dict["rocket_total_mass"].append(rocket_total_mass)
        self.sim_dict["motor_mass"].append(motor_mass)

    # Converts the data saved in this sim into csv
    def to_csv(self):
        #Output
        record = []
        for point in range(len(self.sim_dict["time"])):
            cur_point = []
            cur_point.append(str(self.sim_dict["time"][point]))
            cur_point += list(map(str, self.sim_dict["pos"][point]))
            cur_point += list(map(str, self.sim_dict["vel"][point]))
            cur_point += list(map(str, self.sim_dict["accel"][point]))
            cur_point += list(map(str, self.sim_dict["ang_pos"][point]))
            cur_point += list(map(str, self.sim_dict["ang_vel"][point]))
            cur_point += list(map(str, self.sim_dict["ang_accel"][point]))
            cur_point += map(str, list([self.sim_dict["alpha"][point]]))
            cur_point += map(str, list([self.sim_dict["rocket_total_mass"][point]]))
            cur_point += map(str, list([self.sim_dict["motor_mass"][point]]))
            cur_point += map(str, list([self.sim_dict["flap_ext"][point]]))
            cur_point += map(str, list([self.sensor_dict["baro_alt"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_accel_x"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_accel_y"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_accel_z"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_ang_pos_x"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_ang_pos_y"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_ang_pos_z"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_gyro_x"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_gyro_y"][point]]))
            cur_point += map(str, list([self.sensor_dict["imu_gyro_z"][point]]))
            cur_point += map(str, list([self.sensor_dict["apogee_estimate"][point]]))
            cur_point += map(str, list(self.kalman_dict["x"][point]))
            cur_point += map(str, list(self.kalman_dict["y"][point]))
            cur_point += map(str, list(self.kalman_dict["z"][point]))
            cur_point += map(str, list(self.kalman_dict["cov_x"][point]))
            cur_point += map(str, list(self.kalman_dict["cov_y"][point]))
            cur_point += map(str, list(self.kalman_dict["cov_z"][point]))
            cur_point += map(str, list(self.kalman_dict["rx"][point]))
            cur_point += map(str, list(self.kalman_dict["ry"][point]))
            cur_point += map(str, list(self.kalman_dict["rz"][point]))

            record.append(cur_point)
        return record

        
# if __name__ == '__main__':
    # rocket = Rocket(motor_mass=5,rocket_dry_mass=2)
    # print(rocket.rocket_total_mass)
