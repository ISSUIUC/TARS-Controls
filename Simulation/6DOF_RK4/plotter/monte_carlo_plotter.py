import matplotlib.pyplot as plt
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.properties as prop

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

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def std_dev(data:np.ndarray):
    return np.std(data)

def plotState(state:str, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
    
    # Get the time stamps from the nominal data
    times = nominal_data[:, indices["time"]]
    # print(times)
    
    
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
    return fig, ax

#TODO: FINISH THIS
def plotStateAndEstimate(state:str, estimate:str, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None):
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        
    times = nominal_data[:, indices["time"]]
    
    for point in monte_carlo_data:
        ax.plot(times[:len(point[indices[estimate]])], point[indices[estimate]], '0.8')
        
    # Find the mean of the Monte carlo data and plot it
    mean_estimate = np.mean(monte_carlo_data[:, indices[estimate], :], axis=0)
    ax.plot(times, mean_estimate, 'k--', label="Mean")
    
    # Calculate std dev and plot 
    stddev = np.array([std_dev(monte_carlo_data[:, indices[estimate], i]) for i in range(len(times))])
    ax.plot(times, mean_estimate + 3*stddev, 'tab:orange', label="3*Std Dev")
    ax.plot(times, mean_estimate - 3*stddev, 'tab:orange')
    
    # Plot nominal trajectory
    mean_nominal = np.mean(monte_carlo_data[:, indices[state], :], axis=0)
    ax.plot(times, mean_nominal, label="Mean true trajectory")
    fig.legend()
    fig.suptitle(state)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(state)
    return fig, ax
    
    
def plotRMSE(state:str, estimate:str, sensors:list, nominal_data:np.ndarray, monte_carlo_data:np.ndarray, fig=None, ax=None):
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
    

if __name__ == "__main__":
    folder_path = os.path.join(os.path.dirname(__file__), prop.monte_carlo_output_folder)
    nominal_data = np.load(os.path.join(folder_path, 'nominal.npy'))
    monte_carlo_data = np.load(os.path.join(folder_path, 'sim_data_0.npy'))

    for path in os.listdir(folder_path):
        if path.startswith("sim_data") and not path.endswith("_0.npy"):
            monte_carlo_data = np.dstack((monte_carlo_data, np.array(np.load(os.path.join(folder_path, path)))))
    
    monte_carlo_data = np.transpose(monte_carlo_data, (2, 1, 0))
    print(monte_carlo_data.shape)
    print(nominal_data.shape)
    fig, ax = plotState("x", nominal_data, monte_carlo_data)
    plotState("kf_x", nominal_data, monte_carlo_data)
    plotRMSE("x", "kf_x", ["accel_x", "baro_alt"], nominal_data, monte_carlo_data)
    plotStateAndEstimate("x", "kf_x", nominal_data, monte_carlo_data)
    plt.show()
    