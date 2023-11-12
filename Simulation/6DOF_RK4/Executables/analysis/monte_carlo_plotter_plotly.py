import numpy as np
import pandas as pd
import os
import sys
import plotly.graph_objects as go
from plotly.graph_objs import Mesh3d

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..','..','Configurations')))

import data_loader as dataloader
config = dataloader.config

# Indices of each state variable in the output array
# kf and rkf values are arrays containing position, velocity and acceleration in each axis
indices = {"time":0,
           "pos_x":1,
           "pos_y":2,
           "pos_z":3,
           "vel_x":4,
           "vel_y":5,
           "vel_z":6,
           "accel_x":7,
           "accel_y":8,
           "accel_z":9,
           "ang_pos_x":10,
           "ang_pos_y":11,
           "ang_pos_z":12,
           "ang_vel_x":13,
           "ang_vel_y":14,
           "ang_vel_z":15,
           "ang_accel_x":16,
           "ang_accel_y":17,
           "ang_accel_z":18,
           "alpha":19,
           "rocket_total_mass":20,
           "motor_mass":21,
           "flap_ext":22,
           "baro_alt":23,
           "imu_accel_x":24,
           "imu_accel_y":25,
           "imu_accel_z":26,
           "imu_ang_pos_x":27,
           "imu_ang_pos_y":28,
           "imu_ang_pos_z":29,
           "imu_gyro_x":30,
           "imu_gyro_y":31,
           "imu_gyro_z":32,
           "apogee_estimate":33,
           "kalman_pos_x":34,
           "kalman_vel_x":35,
           "kalman_accel_x":36,
           "kalman_pos_y":37,
           "kalman_vel_y":38,
           "kalman_accel_y":39,
           "kalman_pos_z":40,
           "kalman_vel_z":41,
           "kalman_accel_z":42,
           "pos_cov_x":43,
           "vel_cov_x":44,
           "accel_cov_x":45,
           "pos_cov_y":46,
           "vel_cov_y":47,
           "accel_cov_y":48,
           "pos_cov_z":49,
           "vel_cov_z":50,
           "accel_cov_z":51,
           "kalman_rpos_x":52,
           "kalman_rvel_x":53,
           "kalman_raccel_x":54,
           "kalman_rpos_y":55,
           "kalman_rvel_y":56,
           "kalman_raccel_y":57,
           "kalman_rpos_z":58,
           "kalman_rvel_z":59,
           "kalman_raccel_z":60}

default_output_folder = 'Output'
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def std_dev(data:np.ndarray):
    return np.std(data)

def plotTrajectory3D(states:list, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, color_data:np.ndarray=None, cmap:str=None, ax=None, output_folder=None):
    fig = go.Figure()
    # Make sure the data is 3D
    if len(states) != 3:
        raise ValueError("Can only plot 3D trajectories for 3 states")
    
    # Create arrays for the data
    nominal_states = np.zeros((len(states), len(nominal_data)))
    monte_carlo_states = np.zeros((len(monte_carlo_data), len(states), len(monte_carlo_data[0][0])))
    
    # Generate the data from the nominal and monte carlo data
    for i, state in enumerate(states):
        nominal_states[i,:] = nominal_data[:, indices[state]]
        monte_carlo_states[:,i,:] = monte_carlo_data[:, indices[state], :]
    fig.add_trace(go.Scatter3d(x=nominal_states[0,:], y=nominal_states[1,:], z=nominal_states[2,:], mode='lines', name='Nominal', line=dict(color=velocity, width=4)))
    for i in range(len(monte_carlo_data)):
        fig.add_trace(go.Scatter3d(x=monte_carlo_states[i,0,:], y=monte_carlo_states[i,1,:], z=monte_carlo_states[i,2,:], mode='lines', name=f'Monte Carlo {i}', line=dict(color=velocity, width=4)))
    
    std_devs = np.zeros(int(len(monte_carlo_states[0][0])/10))
    for i in range(0, len(monte_carlo_states[0][0]), 10):
        # Generate points on a circle
        for j in range(len(states)):
            std_devs[int(i/10)] += (3*std_dev(monte_carlo_states[:,j,i]))**2
        std_devs[int(i/10)] = np.sqrt(std_devs[int(i/10)])

    # Iterate over each timestep in the standard deviation matrix
    for i in range(len(std_devs)):
        phi = np.linspace(0, 2*np.pi)
        theta = np.linspace(-np.pi/2, np.pi/2)
        phi, theta=np.meshgrid(phi, theta)

        xs = ((std_devs[0]*np.cos(theta) * np.sin(phi)) + nominal_states[0,i*10])
        ys = (std_devs[1]*np.cos(theta) * np.cos(phi) + nominal_states[1,i*10])
        zs = (std_devs[2]*np.sin(theta) + nominal_states[2,i*10])
        fig.add_surface(x=xs, y=ys, z=zs, showscale=False)

    fig.update_layout(
        title="Monte Python Trajectories",
        scene = dict(
            xaxis_title='X [m]',
            yaxis_title='Y [m]',
            zaxis_title='Z [m]'),
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="RebeccaPurple"
        )
    )
    fig.show()
    
if __name__ == "__main__":
    folder_path = os.path.join(os.path.dirname(__file__), config["meta"]["monte_carlo_output_folder"])
    nominal_data = pd.read_csv(os.path.join(folder_path, 'SimData/nominal.csv')).to_numpy()
    monte_carlo_data = pd.read_csv(os.path.join(folder_path, 'SimData/sim_data_0.csv')).to_numpy()

    for path in os.listdir(os.path.join(folder_path, 'SimData')):
        if path.startswith("sim_data") and not path.endswith("_0.csv"):
            monte_carlo_data = np.dstack((monte_carlo_data, np.array(pd.read_csv(os.path.join(folder_path, 'SimData', path)))))
    
    # Reshape the data to make it easier to work with (Different samples on first axis, different states on second axis, timesteps on third axis)
    monte_carlo_data = np.transpose(monte_carlo_data, (2, 1, 0))

    # Define the output folder as the 'figures' folder in the monte carlo output folder
    output_folder = os.path.join(folder_path, 'figures')
    # fig, ax = plotState("pos_x", nominal_data, monte_carlo_data, output_folder=output_folder)
    # plotState("kalman_pos_x", nominal_data, monte_carlo_data, output_folder=output_folder)
    # plotRMSE("pos_x", "kalman_pos_x", ["accel_x", "baro_alt"], nominal_data, monte_carlo_data, output_folder=output_folder)
    # plotStateAndEstimate("pos_x", "kalman_pos_x", nominal_data, monte_carlo_data, output_folder=output_folder)
    
    # Plot 3D trajectory
    velocity = np.sqrt(nominal_data[:, indices["vel_x"]]**2 + nominal_data[:, indices["vel_y"]]**2 + nominal_data[:, indices["vel_z"]]**2)
    plotTrajectory3D(["pos_z", "pos_y", "pos_x"], nominal_data, monte_carlo_data, velocity, "cool", output_folder=output_folder)