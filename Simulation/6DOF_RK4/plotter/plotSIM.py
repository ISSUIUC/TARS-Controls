import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.properties as prop

def plotter(sim_dict, sensor_dict=0, kalman_dict=0):
    """Plots the simulation data from the simulation dictionary
    
    Args:
        sim_dict (dict): Dictionary containing simulation data
        sensor_dict (dict, optional): Dictionary containing sensor data. Defaults to 0.
        kalmann_dict (dict, optional): Dictionary containing kalman filter data. Defaults to 0.
    """
    fig_linear,(pos_nc,vel_nc,accel_nc) = plt.subplots(3,1,figsize=(15,10), sharex=True)
    fig_linear.suptitle("PYSIM 6DOF LINEAR PLOT", color='#F5B14C', fontsize = 25)

    # Altitude Measurements vs Real Altitude vs Kalman Filter Graph (No Control)
    plt.xlabel("Time (s)", fontsize = 14)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,0], label="X", color="tab:red", linewidth = 2)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,1], label="Y", color="tab:green", linewidth = 2)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,2], label="Z", color="tab:blue", linewidth = 2)
    pos_nc.plot(sensor_dict["time"], sensor_dict["baro_alt"], label="Barometric Altimeter", color="darkred", linestyle = ":",linewidth = 2)
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2)
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2)
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2)
    pos_nc.plot(sensor_dict["time"], sensor_dict["apogee_estimate"], label="Apogee Estimator",color="brown", linewidth = 2)
    pos_nc.axhline(y=sim_dict["pos"][-1,0], label="Simulated Apogee", linestyle = "dashed", color="gray", linewidth = 2)
    pos_nc.axhline(y=prop.des_apogee, label="Desired Apogee", linestyle = "dashed", color="orange", linewidth = 2)
    pos_nc.set_ylabel("Position (m)")
    pos_nc.legend(loc='best')

    # Real Velocity vs Kalman Filter Graph (No Control)
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,0],label="X",color="tab:red", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2)
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2)
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2)
    vel_nc.set_ylabel("Velocity (m/s)")
    vel_nc.legend(loc='best')

    # Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (No Control)
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,0],label="X",color="tab:red", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_x"], label="IMU X", color="darkred", linestyle = ":",linewidth = 2)
    accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_y"], label="IMU Y", color="darkgreen", linestyle = ":",linewidth = 2)
    accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_z"], label="IMU Z", color="darkblue", linestyle = ":",linewidth = 2)
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2)
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2)
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2)
    accel_nc.set_ylabel("Acceleration (m/s^2)")
    accel_nc.legend(loc='best')
    plt.tight_layout()

    fig_angular,(ang_pos_nc, ang_vel_nc, ang_accel_nc, alpha_nc) = plt.subplots(4,1,figsize=(15,10), sharex=True)
    fig_angular.suptitle("PYSIM 6DOF ANGULAR PLOT", color='#F5B14C', fontsize = 25)

    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,0], label="Roll", color="tab:red", linewidth = 2)
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,1], label="Pitch", color="tab:green", linewidth = 2)
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,2], label="Yaw", color="tab:blue", linewidth = 2)
    ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_x"], label="IMU Roll", color="darkred", linestyle = ":",linewidth = 2)
    ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_y"], label="IMU Pitch", color="darkgreen", linestyle = ":",linewidth = 2)
    ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_z"], label="IMU Yaw", color="darkblue", linestyle = ":",linewidth = 2)
    ang_pos_nc.set_ylabel("Position (rad)");  
    ang_pos_nc.legend(loc='best')

    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,0],label="Roll",color="tab:red", linewidth = 2);   
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,2],label="Yaw",color="tab:blue", linewidth = 2);  
    ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_x"], label="IMU Roll Rate", color="darkred", linestyle = ":",linewidth = 2)
    ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_y"], label="IMU Pitch Rate", color="darkgreen", linestyle = ":",linewidth = 2)
    ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_z"], label="IMU Yaw Rate", color="darkblue", linestyle = ":",linewidth = 2)
    ang_vel_nc.set_ylabel("Velocity (rad/s)");  
    ang_vel_nc.legend(loc='best')

    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,0],label="Roll",color="tab:red", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,2],label="Yaw",color="tab:blue", linewidth = 2);    
    ang_accel_nc.set_ylabel("Acceleration (rad\s^2)");   
    ang_accel_nc.legend(loc='best')

    alpha_nc.plot(sim_dict["time"], np.degrees(sim_dict["alpha"]),label="Alpha",color="tab:green", linewidth = 2);    
    alpha_nc.axhline(88 ,label="Alpha Limit",color="tab:red",linewidth = 2);  
    alpha_nc.set_ylabel("Angle of Attack (degrees)");   
    alpha_nc.legend(loc='best')
    plt.tight_layout()

    plt.figure()
    plt.plot(sim_dict["time"], sim_dict["flap_ext"],label="Flap Extension",color="tab:green", linewidth = 2);   
    plt.xlabel("Time (sec)"); 
    plt.ylabel("Flap Extension (m)");  

    plt.show()

if __name__ == "__main__":
    sim_dict = {}
    sensor_dict = {}
    kalman_dict = {}
    output_file = os.path.join(os.path.dirname(__file__), prop.output_file)
    file_data = pd.read_csv(output_file)
    sim_dict["time"] = file_data["time"].values
    sensor_dict["time"] = file_data["time"].values
    kalman_dict["time"] = file_data["time"].values

    for attr in ["pos", "vel", "accel", "ang_pos", "ang_vel", "ang_accel"]:
        sim_dict[attr] = np.array(
            list(zip(file_data[f"{attr}_x"].values, file_data[f"{attr}_y"].values, file_data[f"{attr}_z"].values)))
    for attr in ["kalman_pos", "kalman_vel", "kalman_accel"]:
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
    
    plotter(sim_dict=sim_dict, sensor_dict=sensor_dict, kalman_dict=kalman_dict)
