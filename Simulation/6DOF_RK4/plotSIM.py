from turtle import color
import numpy as np
import matplotlib.pyplot as plt


def plotter(sim_dict, sim_dict_noisy=0, sim_dict_kalman=0, sim_error=0, apogee=0, coast_start=0, sim_runtime=0):

    plt.style.use('Solarize_Light2')


    fig,(alt_nc,vel_nc,accel_nc) = plt.subplots(3,1,figsize=(15,10), sharex=False)
    fig.suptitle("SIM PLOT", color='#F5B14C', fontsize = 30)

    # Altitude Measurements vs Real Altitude vs Kalman Filter Graph (No Control)
    # plt.legend(fontsize = 14) 
    plt.xlabel("Time (s)", fontsize = 14)
    alt_nc.plot(sim_dict["time"], sim_dict["pos"], label="Altitude", color="tab:red", linewidth = 2)
    # alt_nc.axhline(y = apogee, color = "tab:brown", linestyle = "dotted", linewidth = 2.5, label="Apogee")
    # alt_nc.plot(sim_dict_noisy["time"], sim_dict_noisy["altitude"],label="Noisy Altitude Reading",color="lightsteelblue", linewidth = 3, linestyle=":")
    # alt_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["altitude"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # alt_nc.set(ylabel = "Altitude (m)")
    # alt_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # alt_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    # alt_nc.text(0.67, 0.05, "Target Apogee: " + str(round(const.DESIRED_APOGEE,2)) + " m", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    # alt_nc.text(0.67, 0.17, "Apogee: " + str(round(apogee,2)) + " m", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    # alt_nc.text(0.67, 0.29, "Sim Runtime: " + str(round(sim_runtime,2)) + " s", verticalalignment='center', horizontalalignment='left', transform=alt_nc.transAxes, color='black', fontsize=10)
    alt_nc.legend()

    # Real Velocity vs Kalman Filter Graph (No Control)
    vel_nc.plot(sim_dict["time"], sim_dict["vel"],label="Velocity",color="royalblue", linewidth = 2);  
    # vel_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["velocity"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # vel_nc.set(ylabel = "Velocity (m/s)")
    # vel_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # vel_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    vel_nc.legend()

    # Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (No Control)
    accel_nc.plot(sim_dict["time"], sim_dict["accel"],label="Acceleration",color="tab:green", linewidth = 2);  
    # accel_nc.plot(sim_dict_noisy["time"], sim_dict_noisy["acceleration"],label="Noisy Accelerometer Reading",color="lightsteelblue", linewidth = 3, linestyle=":")
    # accel_nc.plot(sim_dict_kalman["time"], sim_dict_kalman["acceleration"],label="Kalman Filter State",linestyle="--",color="tab:red")
    # accel_nc.set(ylabel = "Acceleration (m/s^2)")
    # accel_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # accel_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    accel_nc.legend()

    # error_nc.plot(sim_error["time"], sim_error["error"],label="Error",color="tab:purple", linewidth = 1)
    # error_nc.set(ylabel = "Altitude Error (m)")
    # error_nc.axvline(x = const.sim_start_delay, color = "tab:green", linestyle = "dotted", linewidth = 2.5, label="Boost")
    # error_nc.axvline(x = coast_start, color = "tab:blue", linestyle = "dotted", linewidth = 2.5, label="Coast")
    # error_nc.legend()

    plt.show()
