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
import dynamics.controller as contr

config = dataloader.config

def append_to_array(array, x, time_stamp, baro_alt, accel, bno_ang_pos, gyro, kalman_filter, kalman_filter_r, alpha, apogee_estimation, rocket_total_mass, motor_mass, flap_ext):
    '''Append data to array
    
    Args:
        array (np.array): array to append to
        x (np.array): state vector of rocket
        time_stamp (float): current time stamp
        baro_alt (float): barometric altitude
        accel (np.array): accelerometer data
        bno_ang_pos (np.array): bno055 angular position
        gyro (np.array): gyroscope data
        kalman_filter (np.ndarray): Linear kalman filter data
        kalman_filter_r (np.ndarray): Rotational kalman filter data
        alpha (float): angle of attack
        apogee_estimation (float): current apogee estimate
        rocket_total_mass (float): current rocket mass
        motor_mass (float): current motor mass
        flap_ext (float): current flap extension
        
    Returns:
        (np.array): Array with new datapoint appended
    '''
    
    datapoint = np.array([])
    # Simulation true values
    datapoint = np.append(datapoint, x[0])
    datapoint = np.append(datapoint, x[1])
    datapoint = np.append(datapoint, x[2])
    datapoint = np.append(datapoint, x[3])
    datapoint = np.append(datapoint, x[4])
    datapoint = np.append(datapoint, x[5])
    datapoint = np.append(datapoint, time_stamp)
    datapoint = np.append(datapoint, flap_ext)
    datapoint = np.append(datapoint, alpha)
    datapoint = np.append(datapoint, rocket_total_mass)
    datapoint = np.append(datapoint, motor_mass)
    
    # Sensor values
    datapoint = np.append(datapoint, baro_alt)
    datapoint = np.append(datapoint, accel[0])
    datapoint = np.append(datapoint, accel[1])
    datapoint = np.append(datapoint, accel[2])
    datapoint = np.append(datapoint, bno_ang_pos[0])
    datapoint = np.append(datapoint, bno_ang_pos[1])
    datapoint = np.append(datapoint, bno_ang_pos[2])
    datapoint = np.append(datapoint, gyro[0])
    datapoint = np.append(datapoint, gyro[1])
    datapoint = np.append(datapoint, gyro[2])
    datapoint = np.append(datapoint, apogee_estimation)
    
    # Linear Kalman filter values
    datapoint = np.append(datapoint, kalman_filter[0:3])
    datapoint = np.append(datapoint, kalman_filter[3:6])
    datapoint = np.append(datapoint, kalman_filter[6:9])
    
    # Rotational Kalman filter values
    datapoint = np.append(datapoint, kalman_filter_r[0:3])
    datapoint = np.append(datapoint, kalman_filter_r[3:6])
    datapoint = np.append(datapoint, kalman_filter_r[6:9])

    return np.vstack((array, datapoint)) if len(array) != 0 else np.array([datapoint])

    
def simulator(x0, dt, sample_number, run_folder, target_size:int, nominal:bool=False, wind_direction_variance_mean:float=0, 
              wind_direction_variance_stddev:float=0.01, 
              wind_magnitude_variance_mean:float=0, 
              wind_magnitude_variance_stddev:float=0.5, 
              enable_direction_variance:bool=None, 
              enable_magnitude_variance:bool=None,
              nominal_wind_direction:np.ndarray=np.array([-1,0,0]), 
              nominal_wind_magnitude:float=0):
    '''Runs the simulation with the given x0 and wind variance parameters
    
    Args:
        x0 (np.ndarray): initial state vector
        dt (float): time step between each iteration in simulation
        wind_direction_variance_mean (float, optional): mean of wind direction variance. Defaults to 0.
        wind_direction_variance_stddev (float, optional): standard deviation of wind direction variance. Defaults to 0.01.
        wind_magnitude_variance_mean (float, optional): mean of wind magnitude variance. Defaults to 0.
        wind_magnitude_variance_stddev (float, optional): standard deviation of wind magnitude variance. Defaults to 0.5.
        nominal_wind_direction (np.ndarray, optional): nominal wind direction. Defaults to np.array([-1,0,0]).
        nominal_wind_magnitude (float, optional): nominal wind magnitude. Defaults to 0.
    '''
    t_start = time.time()
    enable_direction_variance = not nominal if enable_direction_variance is None else enable_direction_variance
    enable_magnitude_variance = not nominal if enable_magnitude_variance is None else enable_magnitude_variance
    atm = atmosphere.Atmosphere(wind_direction_variance_mean=wind_direction_variance_mean, 
                                wind_direction_variance_stddev=wind_direction_variance_stddev, 
                                wind_magnitude_variance_mean=wind_magnitude_variance_mean, 
                                wind_magnitude_variance_stddev=wind_magnitude_variance_stddev, 
                                enable_direction_variance=enable_direction_variance, 
                                enable_magnitude_variance=enable_magnitude_variance,
                                nominal_wind_direction=nominal_wind_direction,
                                nominal_wind_magnitude=nominal_wind_magnitude)
    
    # rocket = rocket_model.Rocket(atm=atm)
    rocket = rocket_model.Rocket(config['rocket']['stages'][0], atm=atm)
    motor = rocket.motor
    sim = sim_class.Simulator(atm=atm, rocket=rocket)
    
    # Array to store the data from this simulation. Each row is a datapoint. This replaces the 3 dictionaries used in main.py 
    # for the performance improvement offered by .npy files over .csv files.
    sim_data = np.array([])

    # Get sensor config data
    sensor_config = rocket.stage_config['sensors']
    
    x = x0.copy()
    # Random variations to initial state
    if not nominal:
        x[0,0] *= np.random.uniform(0.9, 1.1)
        x[0,1] *= np.random.uniform(0.9, 1.1)
        x[0,2] *= np.random.uniform(0.9, 1.1)
            
    kalman_filter = ekf.KalmanFilter(
        dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    accel = sensors.get_accelerometer_data(x, sensor_config)
    x_data = []
    y_data = []
    z_data = []
    for i in range(10):
        reading = sensors.get_accelerometer_data(x, sensor_config)
        x_data.append(reading[0])
        y_data.append(reading[1])
        z_data.append(reading[2])
        
    ax = sum(x_data)/len(x_data)
    ay = sum(y_data)/len(y_data)
    az = sum(z_data)/len(z_data)

    pitch = -1 * (np.arctan2(-az,-ay) + np.pi/2)
    yaw = np.arctan2(-ax,-ay) + np.pi/2
    r_kalman_filter = r_ekf.KalmanFilter_R(
        dt, 0.0, 0.0, 0.0, pitch, 0.0, 0.0, yaw, 0.0, 0.0)
    # r_kalman_filter = r_ekf.KalmanFilter_R(
    #     dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    time_stamp = 0

    # Use an n value (last parameter) that is divisible by 3 to make computations easier
    apogee_estimator = apg.Apogee(kalman_filter.get_state(), 0.1, 0.01, 3, 30, atm, rocket.stage_config)
    Kp, Ki, Kd = 0.0002, 0, 0
    controller = contr.Controller(Kp, Ki, Kd, dt, config['desired_apogee'], rocket.stage_config)

    # Idle stage
    while time_stamp < rocket.delay:
        time_stamp += dt
        baro_alt = sensors.get_barometer_data(x, sensor_config)
        accel = sensors.get_accelerometer_data(x, sensor_config)
        gyro = sensors.get_gyro_data(x, sensor_config)
        bno_ang_pos = sensors.get_bno_orientation(x, sensor_config)

        kalman_filter.priori(np.array([0.0, 0.0, 0.0, 0.0]))
        kalman_filter.update(bno_ang_pos, baro_alt,
                             accel[0], accel[1], accel[2])
        
        r_kalman_filter.priori()
        r_kalman_filter.update(*gyro, *accel)

        kalman_filter.reset_lateral_pos()
        current_state = kalman_filter.get_state()
        current_state_r = r_kalman_filter.get_state()
        x = np.round(x, 3)
        time_stamp = np.round(time_stamp, 3)
        baro_alt = np.round(baro_alt, 3)
        accel = np.round(accel, 3)
        bno_ang_pos = np.round(bno_ang_pos, 3)
        gyro = np.round(gyro, 3)
        current_state = np.round(current_state, 3)
        current_state_r = np.round(current_state_r, 3)
        rocket.rocket_total_mass = np.round(rocket.rocket_total_mass, 3)
        rocket.motor_mass = np.round(rocket.motor_mass, 3)
        sim_data = append_to_array(sim_data, x, time_stamp, baro_alt, accel, bno_ang_pos, gyro, current_state, current_state_r, 0, current_state[0], rocket.rocket_total_mass, rocket.motor_mass, 0)


    # print("Ignition")

    # # while x[1][prop.vertical] > prop.apogee_thresh and x[0][prop.vertical] > prop.start_thresh:
    start = True

    while x[1, 0] >= 0 or start:
        if start:
            start = False
        # print("Timestamp: ", time_stamp)
        # Get sensor data
        baro_alt = sensors.get_barometer_data(x, sensor_config)
        accel = sensors.get_accelerometer_data(x, sensor_config)
        gyro = sensors.get_gyro_data(x, sensor_config)
        bno_ang_pos = sensors.get_bno_orientation(x, sensor_config)

        # Kalman Filter stuff goes here
        kalman_filter.priori(np.array([0.0, 0.0, 0.0, 0.0]))
        kalman_filter.update(bno_ang_pos, baro_alt,
                             accel[0], accel[1], accel[2])

        r_kalman_filter.priori()
        r_kalman_filter.update(*gyro, *accel)

        current_state = kalman_filter.get_state()
        current_state_r = r_kalman_filter.get_state()

        apogee_est = apogee_estimator.predict_apogee(current_state[0:3])

        flap_ext = controller.get_flap_extension(time_stamp > rocket.stage_config["motor"]["delay"] and np.linalg.norm(motor.get_thrust(time_stamp)) <= 0, apogee_est)

        # flap_ext will be passed by kalman filter
        rocket.set_motor_mass(time_stamp)
        # rocket.motor_mass = motor.get_mass(time_stamp)

        x, alpha = sim.RK4(x, dt, time_stamp, flap_ext, density_noise=True)
        time_stamp += dt

        x = np.round(x, 3)
        time_stamp = np.round(time_stamp, 3)
        baro_alt = np.round(baro_alt, 3)
        accel = np.round(accel, 3)
        bno_ang_pos = np.round(bno_ang_pos, 3)
        gyro = np.round(gyro, 3)
        current_state = np.round(current_state, 3)
        current_state_r = np.round(current_state_r, 3)
        alpha = np.round(alpha, 3)
        apogee_est = np.round(apogee_est, 3)
        rocket.rocket_total_mass = np.round(rocket.rocket_total_mass, 3)
        rocket.motor_mass = np.round(rocket.motor_mass, 3)
        flap_ext = np.round(flap_ext, 3)
        
        sim_data = append_to_array(sim_data, x, time_stamp, baro_alt, accel, bno_ang_pos, gyro, current_state, current_state_r, alpha, apogee_est, rocket.rocket_total_mass, rocket.motor_mass, flap_ext)

    # Post processing to make sure all datasets are the same size
    if len(sim_data) > target_size:
        sim_data = sim_data[:target_size]
    elif len(sim_data) < target_size:
        while len(sim_data) < target_size:
            sim_data = np.append(sim_data, np.array([sim_data[-1]]), axis=0)
            
    print(sim_data.shape)
    filename = 'nominal.npy' if nominal else f'sim_data_{sample_number}.npy'
    np.save(f"{run_folder}/SimData/{filename}", sim_data)
    t_end = time.time() - t_start
    print(f"Time: {t_end:.2f}")
    
    
def run(x0, dt, samples:int, run_folder:str, target_size:int, clear_contents:bool=False, 
        wind_direction_variance_mean:float=0.0, 
        wind_direction_variance_stddev:float=0.01,
        wind_magnitude_variance_mean:float=0.0,
        wind_magnitude_variance_stddev:float=0.5,
        enable_direction_variance:bool=None,
        enable_magnitude_variance:bool=None,
        nominal_wind_direction:np.ndarray=np.array([-1.0, 0.0, 0.0]),
        nominal_wind_magnitude:float=0.0):
    if clear_contents and os.path.exists(run_folder):
        shutil.rmtree(run_folder)
    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
        os.mkdir(f"{run_folder}/figures")
        os.mkdir(f"{run_folder}/SimData")
    
    # Calculate nominal trajectory
    print("Calculating nominal trajectory")
    simulator(x0, dt, 0, run_folder, target_size, nominal=True, 
              nominal_wind_direction=nominal_wind_direction, 
              nominal_wind_magnitude=nominal_wind_magnitude)
    print("Running Monte Carlo simulations")
    for i in range(samples):
        print(f"Running sample {i+1} of {samples}")
        simulator(x0, dt, i, run_folder, target_size, wind_direction_variance_mean=wind_direction_variance_mean,
                  wind_direction_variance_stddev=wind_direction_variance_stddev,
                  wind_magnitude_variance_mean=wind_magnitude_variance_mean,
                  wind_magnitude_variance_stddev=wind_magnitude_variance_stddev,
                  enable_direction_variance=enable_direction_variance,
                  enable_magnitude_variance=enable_magnitude_variance,
                  nominal_wind_direction=nominal_wind_direction,
                  nominal_wind_magnitude=nominal_wind_magnitude)
        
    print("Done")
    
if __name__ == "__main__":
    '''Main function for running the simulator
    
    Set the initial conditions for the simulation in x0 using the following scheme:
             x:         y:         z:
           [[pos,       pos,       pos],
            [vel,       vel,       vel],
            [accel,     accel,     accel],
            [ang_pos,   ang_pos,   ang_pos],
            [ang_vel,   ang_vel,   ang_vel],
            [ang_accel, ang_accel, ang_accel]]
    
    The time step is set by dt, and the number of samples is set by samples. The timestep is set to 0.1s instead of the default
    standard 0.01s used in main.py to speed up the simulation for testing. For higher fidelity simulations, use 0.01s. Similarly,
    the number of samples is set to 30 instead of something higher (100-1000) to speed up the simulation. Future work will include
    adding multiprocessing to speed up the simulation when using larger sample sizes.
    
    Target size is the number of data points to be saved for each sample. This ensures that each 
    sample has the same number of data points.
    
    The output folder is set by output_folder. This is where the data will be saved. There is a folder
    saved in properties that can be used by default, but you can also change it if you want to save multiple runs with
    different parameters.
    
    The nominal trajectory is defined as the trajectory that the rocket would take if there were no wind. This is useful for quantifying
    trajectory uncertainty, but not as useful for quantifying sensor noise. To analyze sensor noise, disable the direction and magnitude
    variance by setting enable_direction_variance and enable_magnitude_variance to False.
    '''
    # Set up the rocket
    x0 = np.zeros((6,3))
    x0[3] = np.array([0.0, 0.05, 0.01])
    dt = 0.1
    target_size = 880
    samples = 30
    # Calculate nominal trajectory
    output_folder = os.path.join(os.path.dirname(__file__), config["meta"]["monte_carlo_output_folder"])
    run(x0, 
        dt, 
        samples, 
        output_folder, 
        target_size, 
        clear_contents=True, 
        wind_direction_variance_stddev=0.2, 
        wind_magnitude_variance_stddev=3.0,
        enable_direction_variance=True,
        enable_magnitude_variance=True, 
        nominal_wind_magnitude=2,
        nominal_wind_direction=np.array([0,1,0]))