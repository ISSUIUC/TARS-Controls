import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.graph_objs.scatter.marker import Line
from plotly import * 
import numpy as np
import numpy.linalg as la
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.properties as prop
import properties.data_loader as dataloader

config = dataloader.config

def rmse(actual, est):
    return np.sqrt(np.sum(np.square(est - actual))/len(actual))

def plotter(sim_dict, sensor_dict=0, kalman_dict=0):
    """Plots the simulation data from the simulation dictionary
    
    Args:
        sim_dict (dict): Dictionary containing simulation data
        sensor_dict (dict, optional): Dictionary containing sensor data. Defaults to 0.
        kalmann_dict (dict, optional): Dictionary containing kalman filter data. Defaults to 0.
    """
    fig_linear,(pos_nc,vel_nc,accel_nc,flap_nc) = plt.subplots(4,1,figsize=(15,10), sharex=True);   
    #fig_linear.suptitle("PYSIM 6DOF LINEAR PLOT", color='#F5B14C', fontsize = 20); 

    trace_pos_0 = go.Scatter(x = sim_dict["time"], y = sim_dict["pos"][:,0], mode = "lines", name = "X");
    trace_pos_1 = go.Scatter(x = sim_dict["time"], y = sim_dict["pos"][:,1], mode = "lines", name = "Y");
    trace3_pos_2 = go.Scatter(x = sim_dict["time"], y = sim_dict["pos"][:,2], mode = "lines", name = "Z");
    trace4_kalman_pos_0 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,0], mode = "lines", name = "X Estimate", line=dict(dash ='dot'));
    trace4_kalman_pos_1 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,1], mode = "lines", name = "Y Estimate", line=dict(dash ='dot'));
    trace4_kalman_pos_2 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,2], mode = "lines", name = "Z Estimate", line=dict(dash ='dot'));
    trace_apogee_estimate = go.Scatter(x = sensor_dict["time"],y = sensor_dict["apogee_estimate"], mode = "lines", name = "Apogee Estimate", line=dict(dash ='dot'));
    

   # trace_vel_0 = go.Scatter(x = sensor_dict["time"], y = sim_dict["vel"][:,0], mode = "lines", name = "X Estimate");
    trace_vel_0 = go.Scatter(x = sim_dict["time"], y = sim_dict["vel"][:,0], mode = "lines", name = "X");
    trace_vel_1 = go.Scatter(x = sim_dict["time"], y = sim_dict["vel"][:,1], mode = "lines",  name = "Y");
    trace_vel_2 = go.Scatter(x = sim_dict["time"], y = sim_dict["vel"][:,2], mode = "lines",  name = "Z");
    trace_kalman_vel_0 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_vel"][:,0], mode = "lines",  name = "X Estimate", line=dict(dash ='dot'));
    trace_kalman_vel_1 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_vel"][:,1], mode = "lines",  name = "Y Estimate", line=dict(dash ='dot'));
    trace_kalman_vel_2 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_vel"][:,2], mode = "lines",  name = "Z Estimate", line=dict(dash ='dot'));


    trace_accel_0 = go.Scatter(x = sim_dict["time"], y = sim_dict["accel"][:,0], mode = "lines", name = "X");
    trace_accel_1 = go.Scatter(x = sim_dict["time"], y = sim_dict["accel"][:,1], mode = "lines", name = "Y");
    trace_accel_2 = go.Scatter(x = sim_dict["time"], y = sim_dict["accel"][:,2], mode = "lines", name = "Z");
    trace_kalman_accel_0 = go.Scatter(x = kalman_dict["time"],y = kalman_dict["kalman_accel"][:,0], mode = "lines", name = "X Estimate", line=dict(dash ='dot'));
    trace_kalman_accel_1 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_accel"][:,1], mode = "lines", name = "Y Estimate", line=dict(dash ='dot'));
    trace_kalman_accel_2 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_accel"][:,2], mode = "lines", name = "Z Estimate", line=dict(dash ='dot'));
   
    trace_ang_pos_0 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_pos"][:,0], mode = "lines");
    trace_ang_pos_1 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_pos"][:,1], mode = "lines");
    trace_ang_pos_2 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_pos"][:,2], mode = "lines");
    trace_kalman_rpos_0 = go.Scatter(x = kalman_dict["time"],y = kalman_dict["kalman_rpos"][:,0], mode = "lines");
    trace_kalman_rpos_1 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_rpos"][:,1], mode = "lines");
    trace_kalman_rpos_2 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_rpos"][:,2], mode = "lines");
    trace_ang_vel_0 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_vel"][:,0], mode = "lines");
    trace_ang_vel_1 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_vel"][:,1], mode = "lines");
    trace_ang_vel_2 = go.Scatter(x = sim_dict["time"], y = sim_dict["ang_vel"][:,2], mode = "lines");
    trace_kalman_pos_0 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,0], mode = "lines");
    trace_kalman_pos_1 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,1], mode = "lines");
    trace_kalman_pos_2 = go.Scatter(x = kalman_dict["time"], y = kalman_dict["kalman_pos"][:,2], mode = "lines"); 
    
  #  trace8 = go.Line(x = sim_dict["time"], y = sim_dict["vel"][:,0]);
  #  trace9 = go.Line(x = sim_dict["time"], y = sim_dict["vel"][:,1]);

    layout1 = go.Layout(yaxis=dict(title="Position"))
    data1 = [trace_pos_0, trace_pos_1, trace3_pos_2, trace4_kalman_pos_0, trace4_kalman_pos_1, trace4_kalman_pos_2,];
    fig1_linear = go.Figure(layout = layout1, data = data1);
    fig1_linear.show();

    layout2 = go.Layout(yaxis=dict(title="Velocity"))
    data2 = [trace_vel_0, trace_vel_1, trace_vel_2, trace_kalman_vel_0, trace_kalman_vel_1, trace_kalman_vel_2];
    fig2_linear = go.Figure(layout = layout2, data = data2);
    fig2_linear.show();

    layout3 = go.Layout(yaxis=dict(title="Acceleration"))
    data3 = [trace_accel_0, trace_accel_1, trace_accel_2, trace_kalman_accel_0, trace_kalman_accel_1, trace_kalman_accel_2];
    fig3_linear = go.Figure(layout = layout3, data = data3);
    fig3_linear.show();

    #checking if 3D plot works
    x = np.array([0,1])
    y = np.array([1,1])
    z = np.array([1,1])
    df = pd.DataFrame({"x": x, "y":y, "z":z})

    fig = px.line_3d(df, x="x", y="y", z="z")
    fig.show()

    #implementing 3D plot
    dposition = pd.DataFrame({"x": trace_pos_0, "y": trace_pos_1, "z":trace3_pos_2})
    fig = px.line_3d(dposition, x = "x", y = "y", z= "z")
    fig.show()
    


    #big_plot = make_subplots(rows = 3, cols = 1, shared_xaxes= True)
    #big_plot.add_trace((trace_pos_0, trace_pos_1), row = 1, col = 1);
    #big_plot.add_trace(trace_pos_1, row = 2, col = 1);
    #big_plot.add_trace(trace3_pos_2, row = 3, col = 1);
    #big_plot.show();

    # Altitude Measurements vs Real Altitude vs Kalman Filter Graph (No Control)
    plt.xlabel("Time (s)", fontsize = 14);  
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,0], label="X", color="tab:red", linewidth = 2); 
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,1], label="Y", color="tab:green", linewidth = 2);   
    pos_nc.plot(sim_dict["time"], sim_dict["pos"][:,2], label="Z", color="tab:blue", linewidth = 2);    
    # pos_nc.plot(sensor_dict["time"], sensor_dict["baro_alt"], label="Baro Alt", color="darkred", linestyle = ":",linewidth = 2);    
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2);   
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2); 
    pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2);  
    pos_nc.plot(sensor_dict["time"], sensor_dict["apogee_estimate"], label="Apogee Estimate",color="brown", linewidth = 2); 
    pos_nc.axhline(y=sim_dict["pos"][-1,0], label="Sim Apogee", linestyle = "dashed", color="gray", linewidth = 2); 
    pos_nc.axhline(y=config["desired_apogee"], label="Desired Apogee", linestyle = "dashed", color="orange", linewidth = 2); 
    pos_nc.set_ylabel("Position (m)", fontsize = 10);   
    pos_nc.legend(fontsize=10, loc='upper left', ncol=3);   

    # Real Velocity vs Kalman Filter Graph (No Control)
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,0],label="X",color="tab:red", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    vel_nc.plot(sim_dict["time"], sim_dict["vel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2);   
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2); 
    vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2);  
    vel_nc.set_ylabel("Velocity (m/s)", fontsize = 10); 
    vel_nc.legend(fontsize=10, loc='upper left', ncol =2);  

    # Acceleration Measurements vs Real Acceleration vs Kalman Filter Graph (No Control)
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,0],label="X",color="tab:red", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,1],label="Y",color="tab:green", linewidth = 2);  
    accel_nc.plot(sim_dict["time"], sim_dict["accel"][:,2],label="Z",color="tab:blue", linewidth = 2);  
    # accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_x"], label="IMU X", color="darkred", linestyle = ":",linewidth = 2);  
    # accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_y"], label="IMU Y", color="darkgreen", linestyle = ":",linewidth = 2);    
    # accel_nc.plot(sensor_dict["time"], sensor_dict["imu_accel_z"], label="IMU Z", color="darkblue", linestyle = ":",linewidth = 2); 
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,0], label="X Estimate", color="purple", linestyle = "dashed",linewidth = 2);   
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,1], label="Y Estimate", color="lime", linestyle = "dashed",linewidth = 2); 
    accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,2], label="Z Estimate", color="skyblue", linestyle = "dashed",linewidth = 2);  
    accel_nc.set_ylabel("Acceleration (m/s$^2$)", fontsize = 10); 
    accel_nc.legend(fontsize=10, loc='upper left', ncol=2); 

    flap_nc.plot(sim_dict["time"], sim_dict["flap_ext"],label="Flap Extension",color="tab:green", linewidth = 2);   
    #skipped this plot ^^^^
    flap_nc.set_ylabel("Flap Extension (m)", fontsize = 10);  
    flap_nc.legend(fontsize=10, loc='upper left',ncol=1);   

    plt.tight_layout(); 

    fig_angular,(ang_pos_nc, ang_vel_nc, ang_accel_nc, alpha_nc) = plt.subplots(4,1,figsize=(15,10), sharex=True);  
    fig_angular.suptitle("PYSIM 6DOF ANGULAR PLOT", color='#F5B14C', fontsize = 20);    
    plt.xlabel("Time (s)", fontsize = 14);  
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,0], label="Roll", color="tab:red", linewidth = 2);  
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,1], label="Pitch", color="tab:green", linewidth = 2);   
    ang_pos_nc.plot(sim_dict["time"], sim_dict["ang_pos"][:,2], label="Yaw", color="tab:blue", linewidth = 2);  
    # ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_x"], label="IMU Roll", color="darkred", linestyle = ":",linewidth = 2);   
    # ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_y"], label="IMU Pitch", color="darkgreen", linestyle = ":",linewidth = 2);    
    # ang_pos_nc.plot(sensor_dict["time"], sensor_dict["imu_ang_pos_z"], label="IMU Yaw", color="darkblue", linestyle = ":",linewidth = 2);   
    ang_pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_rpos"][:,0], label="Roll Estimate", color="purple", linestyle = "dashed",linewidth = 2);  
    ang_pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_rpos"][:,1], label="Pitch Estimate", color="lime", linestyle = "dashed",linewidth = 2);    
    ang_pos_nc.plot(kalman_dict["time"], kalman_dict["kalman_rpos"][:,2], label="Yaw Estimate", color="skyblue", linestyle = "dashed",linewidth = 2); 
    ang_pos_nc.set_ylabel("Position (rad)", fontsize = 10)
    ang_pos_nc.legend(fontsize=10, loc='upper left',ncol=2)

    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,0],label="Roll",color="tab:red", linewidth = 2);   
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_vel_nc.plot(sim_dict["time"], sim_dict["ang_vel"][:,2],label="Yaw",color="tab:blue", linewidth = 2); 
    # ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_x"], label="IMU Roll Rate", color="darkred", linestyle = ":",linewidth = 2); 
    # ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_y"], label="IMU Pitch Rate", color="darkgreen", linestyle = ":",linewidth = 2);  
    # ang_vel_nc.plot(sensor_dict["time"], sensor_dict["imu_gyro_z"], label="IMU Yaw Rate", color="darkblue", linestyle = ":",linewidth = 2); 
    ang_vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_rvel"][:,0], label="Roll Estimate", color="purple", linestyle = "dashed",linewidth = 2);  
    ang_vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_rvel"][:,1], label="Pitch Estimate", color="lime", linestyle = "dashed",linewidth = 2);    
    ang_vel_nc.plot(kalman_dict["time"], kalman_dict["kalman_rvel"][:,2], label="Yaw Estimate", color="skyblue", linestyle = "dashed",linewidth = 2);  
    ang_vel_nc.set_ylabel("Velocity (rad/s)", fontsize = 10);  
    ang_vel_nc.legend(fontsize=10, loc='upper left',ncol=2)

    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,0],label="Roll",color="tab:red", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,1],label="Pitch",color="tab:green", linewidth = 2);  
    ang_accel_nc.plot(sim_dict["time"], sim_dict["ang_accel"][:,2],label="Yaw",color="tab:blue", linewidth = 2);    
    ang_accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_raccel"][:,0], label="Roll Estimate", color="purple", linestyle = "dashed",linewidth = 2);  
    ang_accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_raccel"][:,1], label="Pitch Estimate", color="lime", linestyle = "dashed",linewidth = 2);    
    ang_accel_nc.plot(kalman_dict["time"], kalman_dict["kalman_raccel"][:,2], label="Yaw Estimate", color="skyblue", linestyle = "dashed",linewidth = 2); 
    ang_accel_nc.set_ylabel("Acceleration (rad/s$^2$)", fontsize = 10);   
    ang_accel_nc.legend(fontsize=10, loc='upper left',ncol=2)

    alpha_nc.plot(sim_dict["time"], np.degrees(sim_dict["alpha"]),label="Alpha",color="tab:green", linewidth = 2);    
    alpha_nc.axhline(88 ,label="Alpha Limit",color="tab:red",linewidth = 2);  
    alpha_nc.set_ylabel("AoA (degrees)", fontsize = 10);   
    alpha_nc.legend(fontsize=10, loc='center left',ncol=1)
    plt.tight_layout() 

    fig_error, (pos_error_nc, vel_error_nc, accel_error_nc) = plt.subplots(3, 1, figsize = (15,10), sharex=True)
    fig_error.suptitle("Kalman Filter Position, Velocity, and Acceleration Error", fontsize = 16)
    plt.xlabel("Time (s)", fontsize = 10)
    pos_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,0] - sim_dict["pos"][:,0], label="X", color="tab:red", linestyle = "solid",linewidth = 2);
    pos_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["pos_cov"][:,0]), label="$\pm3\sigma$", color="tab:red", linestyle = "dashed",linewidth = 1);
    pos_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["pos_cov"][:,0]), color="tab:red", linestyle = "dashed",linewidth = 1);
    pos_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,1] - sim_dict["pos"][:,1], label="Y", color="tab:green", linestyle = "solid",linewidth = 2);
    pos_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["pos_cov"][:,1]), label="$\pm3\sigma$", color="tab:green", linestyle = "dashed",linewidth = 1);
    pos_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["pos_cov"][:,1]), color="tab:green", linestyle = "dashed",linewidth = 1);
    pos_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_pos"][:,2] - sim_dict["pos"][:,2], label="Z", color="tab:blue", linestyle = "solid",linewidth = 2);
    pos_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["pos_cov"][:,2]), label="$\pm3\sigma$", color="tab:blue", linestyle = "dashed",linewidth = 1);
    pos_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["pos_cov"][:,2]), color="tab:blue", linestyle = "dashed",linewidth = 1);
    pos_error_nc.set_ylabel("Position Error (m)", fontsize = 10)
    pos_error_nc.legend(fontsize=10, loc='upper left',ncol=3)

    vel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,0] - sim_dict["vel"][:,0], label="X", color="tab:red", linestyle = "solid",linewidth = 2);
    vel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["vel_cov"][:,0]), label="$\pm3\sigma$", color="tab:red", linestyle = "dashed",linewidth = 1);
    vel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["vel_cov"][:,0]), color="tab:red", linestyle = "dashed",linewidth = 1);
    vel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,1] - sim_dict["vel"][:,1], label="Y", color="tab:green", linestyle = "solid",linewidth = 2);
    vel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["vel_cov"][:,1]), label="$\pm3\sigma$", color="tab:green", linestyle = "dashed",linewidth = 1);
    vel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["vel_cov"][:,1]), color="tab:green", linestyle = "dashed",linewidth = 1);
    vel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_vel"][:,2] - sim_dict["vel"][:,2], label="Z", color="tab:blue", linestyle = "solid",linewidth = 2);
    vel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["vel_cov"][:,2]), label="$\pm3\sigma$", color="tab:blue", linestyle = "dashed",linewidth = 1);
    vel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["vel_cov"][:,2]), color="tab:blue", linestyle = "dashed",linewidth = 1);
    vel_error_nc.set_ylabel("Velocity Error (m/s)", fontsize = 10)
    vel_error_nc.legend(fontsize=10, loc='upper left',ncol=3)

    accel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,0] - sim_dict["accel"][:,0], label="X", color="tab:red", linestyle = "solid",linewidth = 2);
    accel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["accel_cov"][:,0]), label="$\pm3\sigma$", color="tab:red", linestyle = "dashed",linewidth = 1);
    accel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["accel_cov"][:,0]), color="tab:red", linestyle = "dashed",linewidth = 1);
    accel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,1] - sim_dict["accel"][:,1], label="Y", color="tab:green", linestyle = "solid",linewidth = 2);
    accel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["accel_cov"][:,1]), label="$\pm3\sigma$", color="tab:green", linestyle = "dashed",linewidth = 1);
    accel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["accel_cov"][:,1]), color="tab:green", linestyle = "dashed",linewidth = 1);
    accel_error_nc.plot(kalman_dict["time"], kalman_dict["kalman_accel"][:,2] - sim_dict["accel"][:,2], label="Z", color="tab:blue", linestyle = "solid",linewidth = 2);
    accel_error_nc.plot(kalman_dict["time"], 3*np.sqrt(kalman_dict["accel_cov"][:,2]), label="$\pm3\sigma$", color="tab:blue", linestyle = "dashed",linewidth = 1);
    accel_error_nc.plot(kalman_dict["time"], -3*np.sqrt(kalman_dict["accel_cov"][:,2]), color="tab:blue", linestyle = "dashed",linewidth = 1);
    accel_error_nc.set_ylabel("Acceleration Error (m/s$^2$)", fontsize = 10)
    accel_error_nc.legend(fontsize=10, loc='upper left',ncol=3)
    
    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    fig_3d.suptitle("3D plots", fontsize=20)
    x, y, z = sim_dict["pos"][:,0], sim_dict["pos"][:,1], sim_dict["pos"][:,2]
    vx, vy, vz = sim_dict["vel"][:,0], sim_dict["vel"][:,1], sim_dict["vel"][:,2]
    velocity_magnitude = np.sqrt(vx**2 + vy**2 + vz**2)
    plot = ax_3d.scatter(y, z, x, label="Simulated Position", c=velocity_magnitude, cmap='viridis', linewidth=0.5)
    ax_3d.set_zlabel("Altitude (m)")
    fig_3d.colorbar(plot)
   # plt.show()

if __name__ == "__main__":
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

    plotter(sim_dict=sim_dict, sensor_dict=sensor_dict, kalman_dict=kalman_dict)