import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..','..','Configurations')))

from data_loader import config

def load_record(file_path):
    sim_dict = {}
    sensor_dict = {}
    kalman_dict = {}
    output_file = os.path.join(os.path.dirname(__file__), config["meta"]["output_file"])
    file_data = pd.read_csv(output_file)
    sim_dict["time"] = file_data["time"].values
    sensor_dict["time"] = file_data["time"].values
    kalman_dict["time"] = file_data["time"].values

    for attr in ["pos", "vel", "accel", "ang_pos", "ang_vel", "ang_accel"]:
        sim_dict[attr] = np.array(
            list(zip(file_data[f"{attr}_x"].values, file_data[f"{attr}_y"].values, file_data[f"{attr}_z"].values)))
    for attr in ["kalman_pos", "kalman_vel", "kalman_accel", "pos_cov", "vel_cov", "accel_cov"]:
        kalman_dict[attr] = np.array(
            list(zip(file_data[f"{attr}_x"].values, file_data[f"{attr}_y"].values, file_data[f"{attr}_z"].values)))
    for attr in ["kalman_rpos", "kalman_rvel", "kalman_raccel"]:
        kalman_dict[attr] = np.array(
            list(zip(file_data[f"{attr}_x"].values, file_data[f"{attr}_y"].values, file_data[f"{attr}_z"].values)))
    sim_dict["alpha"] = np.array(list(zip(file_data["alpha"].values)))
    sim_dict["rocket_total_mass"] = np.array(list(zip(file_data["rocket_total_mass"].values)))
    sim_dict["motor_mass"] = np.array(list(zip(file_data["motor_mass"].values)))
    sim_dict["flap_ext"] = np.array(list(zip(file_data["flap_ext"].values)))
    sensor_dict["baro_alt"] = np.array(list(zip(file_data["baro_alt"].values)))
    sensor_dict["imu_accel_x"] = np.array(list(zip(file_data["imu_accel_x"].values)))
    sensor_dict["imu_accel_y"] = np.array(list(zip(file_data["imu_accel_y"].values)))
    sensor_dict["imu_accel_z"] = np.array(list(zip(file_data["imu_accel_z"].values)))
    sensor_dict["imu_ang_pos_x"] = np.array(list(zip(file_data["imu_ang_pos_x"].values)))
    sensor_dict["imu_ang_pos_y"] = np.array(list(zip(file_data["imu_ang_pos_y"].values)))
    sensor_dict["imu_ang_pos_z"] = np.array(list(zip(file_data["imu_ang_pos_z"].values)))
    sensor_dict["imu_gyro_x"] = np.array(list(zip(file_data["imu_gyro_x"].values)))
    sensor_dict["imu_gyro_y"] = np.array(list(zip(file_data["imu_gyro_y"].values)))
    sensor_dict["imu_gyro_z"] = np.array(list(zip(file_data["imu_gyro_z"].values)))
    sensor_dict["apogee_estimate"] = np.array(list(zip(file_data["apogee_estimate"].values)))

    return sim_dict, sensor_dict, kalman_dict

def save_record(record, output_file, target_size=None):
    print("Writing to file...")
    if target_size is not None:
        if len(record) > target_size:
                record = record[:target_size]
        elif len(record) < target_size:
            while len(record) < target_size:
                record = np.append(record, np.array([record[-1]]), axis=0)
    
    with open(output_file, 'w') as f:
        f.write("time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,accel_x,accel_y,accel_z,ang_pos_x,ang_pos_y,ang_pos_z,ang_vel_x,ang_vel_y,ang_vel_z,ang_accel_x,ang_accel_y,ang_accel_z,alpha,rocket_total_mass,motor_mass,flap_ext,baro_alt,imu_accel_x,imu_accel_y,imu_accel_z,imu_ang_pos_x,imu_ang_pos_y,imu_ang_pos_z,imu_gyro_x,imu_gyro_y,imu_gyro_z,apogee_estimate,kalman_pos_x,kalman_vel_x,kalman_accel_x,kalman_pos_y,kalman_vel_y,kalman_accel_y,kalman_pos_z,kalman_vel_z,kalman_accel_z,pos_cov_x,vel_cov_x,accel_cov_x,pos_cov_y,vel_cov_y,accel_cov_y,pos_cov_z,vel_cov_z,accel_cov_z,kalman_rpos_x,kalman_rvel_x,kalman_raccel_x,kalman_rpos_y,kalman_rvel_y,kalman_raccel_y,kalman_rpos_z,kalman_rvel_z,kalman_raccel_z\n")
        for point in record:
            f.write(f"{','.join(point)}\n")
    