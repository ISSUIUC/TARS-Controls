import matplotlib.pyplot as plt
import numpy as np
import os
import sys

plt.rcParams.update({'font.size': 14})

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.properties as prop
import properties.data_loader as dataloader
config = dataloader.config

# Indices of each state variable in the output array
# kf and rkf values are arrays containing position, velocity and acceleration in each axis
indices = {"x":0,
           "y":1,
           "z":2,
           "vx":3,
           "vy":4,
           "vz":5,
           "ax":6,
           "ay":7,
           "az":8,
           "roll":9,
           "pitch":10,
           "yaw":11,
           "ang_vel_roll":12,
           "ang_vel_pitch":13,
           "ang_vel_yaw":14,
           "ang_accel_roll":15,
           "ang_accel_pitch":16,
           "ang_accel_yaw":17,
           "time":18,
           "flap_ext":19,
           "alpha":20,
           "rocket_total_mass":21,
           "motor_mass":22,
           "baro_alt":23,
           "accel_x":24,
           "accel_y":25,
           "accel_z":26,
           "imu_x":27,
           "imu_y":28,
           "imu_z":29,
           "gyro_x":30,
           "gyro_y":31,
           "gyro_z":32,
           "apogee_estimate":33,
           "kf_x":34,
           "kf_vx":35,
           "kf_ax":36,
           "kf_y":37,
           "kf_vy":38,
           "kf_ay":39,
           "kf_z":40,
           "kf_vz":41,
           "kf_az":42,
           "rkf_r":43,
           "rkf_vr":44,
           "rkf_ar":45,
           "rkf_p":46,
           "rkf_vp":47,
           "rkf_ap":48,
           "rkf_y":49,
           "rkf_vy":50,
           "rkf_ay":51}

default_output_folder = 'Output'
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def std_dev(data:np.ndarray):
    return np.std(data)

def plotState(state:str, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None, output_folder=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
    
    # Get the time stamps from the nominal data
    times = nominal_data[:, indices["time"]]    
    
    # Plot the Monte Carlo data
    for point in monte_carlo_data:
        ax.plot(times[:len(point[indices[state]])], point[indices[state]], '0.8')
    
    # Find the mean of the Monte carlo data and plot it
    mean = np.mean(monte_carlo_data[:, indices[state], :], axis=0)
    ax.plot(times, mean, 'k--', label="Mean")
    
    # Calculate std dev and plot 
    stddev = np.array([std_dev(monte_carlo_data[:, indices[state], i]) for i in range(len(times))])
    ax.plot(times, mean + 3*stddev, 'tab:orange', label="3*Std Dev")
    ax.plot(times, mean - 3*stddev, 'tab:orange')
    
    # Plot the nominal data
    ax.plot(times, nominal_data[:, indices[state]], label="Nominal")
    fig.legend()
    fig.suptitle(state)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(state)
    if output_folder is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), '../', output_folder, state + ".png"))
    return fig, ax

#TODO: FINISH THIS
def plotStateAndEstimate(state:str, estimate:str, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None, output_folder=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        
    times = nominal_data[:, indices["time"]]
    
    for point in monte_carlo_data:
        ax.plot(times[:len(point[indices[estimate]])], point[indices[estimate]], '0.8')
        
    # Find the mean of the Monte carlo data and plot it
    mean_estimate = np.mean(monte_carlo_data[:, indices[estimate], :], axis=0)
    ax.plot(times, mean_estimate, 'k--', label="Mean Monte Carlo")
    
    # Calculate std dev and plot 
    stddev = np.array([std_dev(monte_carlo_data[:, indices[estimate], i]) for i in range(len(times))])
    ax.plot(times, mean_estimate + 3*stddev, 'tab:orange', label="3*Std Dev")
    ax.plot(times, mean_estimate - 3*stddev, 'tab:orange')
    
    # Plot nominal trajectory
    mean_nominal = np.mean(monte_carlo_data[:, indices[state], :], axis=0)
    ax.plot(times, mean_nominal, label="Mean true trajectory")
    fig.legend()
    fig.suptitle(state + ' true and estimate')
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(state)
    if output_folder is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), '../', output_folder, state + "_estimate.png"))
    return fig, ax
    
    
def plotRMSE(state:str, estimate:str, sensors:list, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None, output_folder=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
    
    # Get the time stamps from the nominal data
    times = nominal_data[:, indices["time"]]
    
    # calculate the RMSE at each time step for the sensors and the estimate
    sensor_rmse = np.zeros(len(times))
    for sensor in sensors:
        sensor_rmse += np.array([rmse(monte_carlo_data[:, indices[sensor], i], nominal_data[i, indices[sensor]]) for i in range(len(times))])
    
    # calculate RMSE of the estimate at each timestep
    estimate_rmse = np.array([rmse(monte_carlo_data[:, indices[estimate], i], nominal_data[i, indices[estimate]]) for i in range(len(times))])
    
    ax.plot(times, sensor_rmse, label="Sensor RMSE")
    ax.plot(times, estimate_rmse, label="Estimate RMSE")
    fig.legend()
    fig.suptitle("RMSE of " + state + " and " + estimate)
    plt.xlabel("Time (s)")
    plt.ylabel("RMSE")
    if output_folder is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), '../', output_folder, state + "_RMSE.png"))
    return fig, ax
    
def plotTrajectory3D(states:list, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, color_data:str=None, cmap:str=None, fig=None, ax=None, output_folder=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel(states[0])
        ax.set_ylabel(states[1])
        ax.set_zlabel(states[2])
    
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
    
    # Plot nominal data
    cbar = fig.colorbar(ax.scatter(nominal_states[0,:], nominal_states[1,:], nominal_states[2,:], c=color_data, cmap=cmap, label="Nominal", marker='.', linewidth=1, location='bottom'), pad=0.1)
    cbar.set_label('Velocity (m/s)')
    # Generate uncertainty cone based on monte carlo data
    # Every 10th point, generate points on a circle with radius 3 times the standard deviation of the monte carlo data that is orthogonal to the nominal trajectory
    # Calculate standard deviation of the monte carlo data at each time step
    std_devs = np.zeros(int(len(monte_carlo_states[0][0])/10))
    for i in range(0, len(monte_carlo_states[0][0]), 10):
        # Generate points on a circle
        for j in range(len(states)):
            std_devs[int(i/10)] += (3*std_dev(monte_carlo_states[:,j,i]))**2
        std_devs[int(i/10)] = np.sqrt(std_devs[int(i/10)])
    
    # Generate unit vectors that span the circle in the plane orthogonal to the nominal trajectory
    angles = np.linspace(0, 2*np.pi, 10)
    
    # Iterate over each timestep in the standard deviation matrix
    for i in range(len(std_devs)):
        # Calculate the unit vector along the nominal trajectory at this timestep
        nominal_unit_vector = np.zeros(len(states))
        for j in range(len(states)):
            nominal_unit_vector[j] = nominal_states[j, 10*i+1] - nominal_states[j, 10*i]
        nominal_unit_vector = nominal_unit_vector / np.linalg.norm(nominal_unit_vector)

        # Calculate the unit vector orthogonal to the nominal trajectory at each angle in angles
        # Orthogonal vector found by setting the x and y coordinates to the cos and sin of the angle and solving for the
        # third coordinate by setting the dot product of the nominal unit vector and the orthogonal unit vector to 0 then normalizing
        circle_points = np.zeros((len(states), len(angles)))
        for j, angle in enumerate(angles):
            orthogonal_unit_vector = np.zeros(len(states))
            orthogonal_unit_vector[0] = np.cos(angle)
            orthogonal_unit_vector[1] = np.sin(angle)
            orthogonal_unit_vector[2] = -(nominal_unit_vector[0]*orthogonal_unit_vector[0] + nominal_unit_vector[1]*orthogonal_unit_vector[1]) / nominal_unit_vector[2]
            orthogonal_unit_vector = orthogonal_unit_vector / np.linalg.norm(orthogonal_unit_vector)
            circle_points[:,j] = nominal_states[:,10*i] + std_devs[i] * orthogonal_unit_vector
        ax.plot3D(circle_points[0,:], circle_points[1,:], circle_points[2,:], 'gray', alpha=0.5)
        
    fig.legend()
    if output_folder is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), '../', output_folder, "_".join(states) + ".png"))
        
if __name__ == "__main__":
    folder_path = os.path.join(os.path.dirname(__file__), config["meta"]["monte_carlo_output_folder"])
    nominal_data = np.load(os.path.join(folder_path, 'SimData/nominal.npy'))
    monte_carlo_data = np.load(os.path.join(folder_path, 'SimData/sim_data_0.npy'))

    for path in os.listdir(os.path.join(folder_path, 'SimData')):
        if path.startswith("sim_data") and not path.endswith("_0.npy"):
            monte_carlo_data = np.dstack((monte_carlo_data, np.array(np.load(os.path.join(folder_path, 'SimData', path)))))
    
    # Reshape the data to make it easier to work with (Different samples on first axis, different states on second axis, timesteps on third axis)
    monte_carlo_data = np.transpose(monte_carlo_data, (2, 1, 0))

    # Define the output folder as the 'figures' folder in the monte carlo output folder
    output_folder = os.path.join(folder_path, 'figures')
    fig, ax = plotState("x", nominal_data, monte_carlo_data, output_folder=output_folder)
    plotState("kf_x", nominal_data, monte_carlo_data, output_folder=output_folder)
    plotRMSE("x", "kf_x", ["accel_x", "baro_alt"], nominal_data, monte_carlo_data, output_folder=output_folder)
    plotStateAndEstimate("x", "kf_x", nominal_data, monte_carlo_data, output_folder=output_folder)
    
    # Plot 3D trajectory
    velocity = np.sqrt(nominal_data[:, indices["vx"]]**2 + nominal_data[:, indices["vy"]]**2 + nominal_data[:, indices["vz"]]**2)
    plotTrajectory3D(["z", "y", "x"], nominal_data, monte_carlo_data, velocity, "cool", output_folder=output_folder)
    plt.show()
    