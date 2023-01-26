import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

import properties as prop

def plotter(sim_dict, sim_dict_noisy=0, sim_dict_kalman=0, sim_error=0, apogee=0, coast_start=0, sim_runtime=0):

    # plt.style.use('Solarize_Light2')


    fig_linear,(pos_nc,vel_nc,accel_nc) = plt.subplots(3,1,figsize=(15,10), sharex=False)
    fig_linear.suptitle("PYSIM 6DOF LINEAR PLOT", color='#F5B14C', fontsize = 25)

    # Altitude Measurements vs Real Altitude vs Kalman Filter Graph (No Control)
    # plt.legend(fontsize = 14) 
    plt.xlabel("Time (s)", fontsize = 14)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,0], label="X", color="tab:red", linewidth = 2)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,1], label="Y", color="tab:green", linewidth = 2)
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,2], label="Z", color="tab:blue", linewidth = 2)
    pos_nc.set_ylabel("Position (m)")
    # alt_nc.axhline(y = apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Apogee")
    # alt_nc.plot(sim_dict_noisy["time"], sim_dict_noisy["altitude"],label="Noisy Altitude Reading",color="lightsteelblue", linewidth = 3, linestyle=":")
    # alt_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["altitude"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # alt_nc.set(ylabel = "Altitude (m)")
    # alt_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # alt_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    # alt_nc.text(0.67, 0.05, "Target Apogee: " + str(round(const.DESIRED_APOGEE,2)) + " m", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    # alt_nc.text(0.67, 0.17, "Apogee: " + str(round(apogee,2)) + " m", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    # alt_nc.text(0.67, 0.29, "Sim Runtime: " + str(round(sim_runtime,2)) + " s", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    pos_nc.legend()

    # Real Velocity vs Kalman Filter Graph (No Control)
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,0],label="X",color="tab:red", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    vel_nc.set_ylabel("Velocity (m/s)")
    # vel_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["velocity"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # vel_nc.set(ylabel = "Velocity (m/s)")
    # vel_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # vel_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    vel_nc.legend()

    # Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (No Control)
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,0],label="X",color="tab:red", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    accel_nc.set_ylabel("Acceleration (m/s^2)")
    # accel_nc.plot(sim_dict_noisy["time"], sim_dict_noisy["acceleration"],label="Noisy Accelerometer Reading",color="lightsteelblue", linewidth = 3, linestyle=":")
    # accel_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["acceleration"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # accel_nc.set(ylabel = "Acceleration (m/s^2)")
    # accel_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # accel_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    accel_nc.legend()
    plt.tight_layout()

    # error_nc.plot(sim_error["time"], sim_error["error"],label="Error",color="tab:purple", linewidth = 1)
    # error_nc.set(ylabel = "Altitude Error (m)")
    # error_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # error_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    # error_nc.legend()

    fig_angular,(ang_pos_nc, ang_vel_nc, ang_accel_nc, alpha_nc) = plt.subplots(4,1,figsize=(15,10), sharex=False)
    fig_angular.suptitle("PYSIM 6DOF ANGULAR PLOT", color='#F5B14C', fontsize = 25)

    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,0], label="Roll", color="tab:red", linewidth = 2)
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,1], label="Pitch", color="tab:green", linewidth = 2)
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,2], label="Yaw", color="tab:blue", linewidth = 2)
    ang_pos_nc.set_ylabel("Position (rad)");  
    ang_pos_nc.legend()

    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,0],label="Roll",color="tab:red", linewidth = 2);   
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,2],label="Yaw",color="tab:blue", linewidth = 2);  
    ang_vel_nc.set_ylabel("Velocity (rad/s)");  
    ang_vel_nc.legend()

    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,0],label="Roll",color="tab:red", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,2],label="Yaw",color="tab:blue", linewidth = 2);    
    ang_accel_nc.set_ylabel("Acceleration (rad\s^2)");   
    ang_accel_nc.legend()

    alpha_nc.plot(sim_dict["time"], np.degrees(sim_dict["alpha"]),label="Alpha",color="tab:green", linewidth = 2);    
    alpha_nc.axhline(15 ,label="Alpha Limit",color="tab:red",linewidth = 2);  
    alpha_nc.set_ylabel("Angle of Attack (degrees)");   
    alpha_nc.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sim_dict = {}
    output_file = os.path.join(os.path.dirname(__file__), prop.output_file)
    file_data = pd.read_csv(output_file)
    sim_dict["time"] = file_data["time"].values
    # print(sim_dict["time"].values)
    # sim_dict["pos"] = np.array(
    #     list(zip(file_data["pos_x"].values, file_data["pos_y"].values, file_data["pos_z"].values)))
    # print(sim_dict["pos"])
    # sim_dict["vel"] = np.array(
    #     list(zip(file_data["vel_x"].values, file_data["vel_y"].values, file_data["vel_z"].values)))
    # sim_dict["accel"] = np.array(
    #     list(zip(file_data["accel_x"].values, file_data["accel_y"].values, file_data["accel_z"].values)))

    for attr in ["pos", "vel", "accel", "ang_pos", "ang_vel", "ang_accel"]:
        sim_dict[attr] = np.array(
            list(zip(file_data[f"{attr}_x"].values, file_data[f"{attr}_y"].values, file_data[f"{attr}_z"].values)))
    sim_dict["alpha"] = np.array(list(zip(file_data["alpha"].values)))
    plotter(sim_dict=sim_dict)
