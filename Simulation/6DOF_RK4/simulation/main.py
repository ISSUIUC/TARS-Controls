# ______      _____ _              __________ ___________ 
# | ___ \    /  ___(_)            / ___|  _  \  _  |  ___|
# | |_/ /   _\ `--. _ _ __ ___   / /___| | | | | | | |_   
# |  __/ | | |`--. \ | '_ ` _ \  | ___ \ | | | | | |  _|  
# | |  | |_| /\__/ / | | | | | | | \_/ | |/ /\ \_/ / |    
# \_|   \__, \____/|_|_| |_| |_| \_____/___/  \___/\_|    
#        __/ |                                            
#       |___/                                             

# A 6DOF RK-4 Based simulation that uses RASAero aerodynamic data and known motor thrust data to simulate motion of the rocket
# with simulated sensor data as well as a implementation of the Extended Kalman Filter and active drag PID controller for the ISS
# Spaceshot entry for the 2023 IREC competition. This simulation was used to quantify the effects of the airbrakes, test 
# different system design methodlogies, and provide preliminary tuning for the EKF and controller prior to implementation in 
# SILSIM and flight software.
# 
# 2022-2023 Guidance, Navigation, and Control Main Contributors #
# Sub-Team Lead: Parth Shrotri (2024)
# Colin Kinsey (2024)
# Evan Yu (2025)
# Rithvik Bhogavilli (2025)
# Kabir Cheema (2025)
# Freya Bansal (2025)
# Ishaan Bansal (2025)
# Ethan Pereira (2026)

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import estimation.ekf as ekf
import estimation.r_ekf as r_ekf
import properties.properties as prop
import simulator as sim_class
import plotter.plotSIM as plotter
import dynamics.sensors as sensors
import time
import estimation.apogee_estimator as apg
import dynamics.rocket as rocket_model
import environment.atmosphere as atmosphere
import dynamics.controller as contr

atm = atmosphere.Atmosphere()
rocket = rocket_model.Rocket(atm=atm)
motor = rocket.motor
sim = sim_class.Simulator(atm=atm, rocket=rocket)
sim_dict = {
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

kalman_dict = {
    "x": [],
    "y": [],
    "z": [],
    "rx": [],
    "ry": [],
    "rz": [],
    "time": []
}

sensor_dict = {
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

def addToDict(x, baro_alt, accel, bno_ang_pos, gyro, kalman_filter, kalman_filter_r, alpha, apogee_estimation, rocket_total_mass, motor_mass, flap_ext):
    # Append to sensor_dict
    sensor_dict["baro_alt"].append(baro_alt)
    sensor_dict["imu_accel_x"].append(accel[0])
    sensor_dict["imu_accel_y"].append(accel[1])
    sensor_dict["imu_accel_z"].append(accel[2])
    sensor_dict["imu_ang_pos_x"].append(bno_ang_pos[0])
    sensor_dict["imu_ang_pos_y"].append(bno_ang_pos[1])
    sensor_dict["imu_ang_pos_z"].append(bno_ang_pos[2])
    sensor_dict["imu_gyro_x"].append(gyro[0])
    sensor_dict["imu_gyro_y"].append(gyro[1])
    sensor_dict["imu_gyro_z"].append(gyro[2])
    sensor_dict["apogee_estimate"].append(apogee_estimation)

    kalman_dict["x"].append(kalman_filter[0:3])
    kalman_dict["y"].append(kalman_filter[3:6])
    kalman_dict["z"].append(kalman_filter[6:9])

    kalman_dict["rx"].append(kalman_filter_r[0:3])
    kalman_dict["ry"].append(kalman_filter_r[3:6])
    kalman_dict["rz"].append(kalman_filter_r[6:9])

    # Update Simulator Log
    sim_dict["pos"].append(x[0])
    sim_dict["vel"].append(x[1])
    sim_dict["accel"].append(x[2])
    sim_dict["ang_pos"].append(x[3])
    sim_dict["ang_vel"].append(x[4])
    sim_dict["ang_accel"].append(x[5])
    sim_dict["time"].append(sim_dict["time"][-1] +
                            dt if len(sim_dict["time"]) > 0 else 0)
    sim_dict["flap_ext"].append(flap_ext)                    
    sim_dict["alpha"].append(alpha)
    sim_dict["rocket_total_mass"].append(rocket_total_mass)
    sim_dict["motor_mass"].append(motor_mass)


def simulator(x0, dt) -> None:
    '''Method which handles running the simulation and logging sim data to dict

    Args:
        x0 (np.array): state vector initialized to 0s [6x3]
             x:         y:         z:
           [[pos,       pos,       pos],
            [vel,       vel,       vel],
            [accel,     accel,     accel],
            [ang_pos,   ang_pos,   ang_pos],
            [ang_vel,   ang_vel,   ang_vel],
            [ang_accel, ang_accel, ang_accel]]
        dt (float): time step between each iteration in simulation

    '''

    x = x0.copy()
    kalman_filter = ekf.KalmanFilter(
        dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    r_kalman_filter = r_ekf.KalmanFilter_R(
        dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    time_stamp = 0

    # Use an n value (last parameter) that is divisible by 3 to make computations easier
    apogee_estimator = apg.Apogee(kalman_filter.get_state(), 0.1, 0.01, 3, 30, atm)
    Kp, Ki, Kd = 0.0002, 0, 0
    controller = contr.Controller(Kp, Ki, Kd, dt, prop.des_apogee)

    # Idle stage
    while time_stamp < rocket.delay:
        time_stamp += dt
        baro_alt = sensors.get_barometer_data(x)
        accel = sensors.get_accelerometer_data(x)
        gyro = sensors.get_gyro_data(x)
        bno_ang_pos = sensors.get_bno_orientation(x)

        kalman_filter.priori(np.array([0.0, 0.0, 0.0, 0.0]))
        kalman_filter.update(bno_ang_pos, baro_alt,
                             accel[0], accel[1], accel[2])
        
        r_kalman_filter.priori()
        r_kalman_filter.update(*gyro, *accel)

        kalman_filter.reset_lateral_pos()
        current_state = kalman_filter.get_state()
        current_state_r = r_kalman_filter.get_state()

        addToDict(x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_state_r, 0, current_state[0], rocket.rocket_total_mass, rocket.motor_mass, 0)


    print("Ignition")

    # # while x[1][prop.vertical] > prop.apogee_thresh and x[0][prop.vertical] > prop.start_thresh:
    start = True
    t_start = time.time()

    while x[1, 0] >= 0 or start:
        if start:
            start = False
        # print("Timestamp: ", time_stamp)
        # Get sensor data
        baro_alt = sensors.get_barometer_data(x)
        accel = sensors.get_accelerometer_data(x)
        gyro = sensors.get_gyro_data(x)
        bno_ang_pos = sensors.get_bno_orientation(x)

        # Kalman Filter stuff goes here
        kalman_filter.priori(np.array([0.0, 0.0, 0.0, 0.0]))
        kalman_filter.update(bno_ang_pos, baro_alt,
                             accel[0], accel[1], accel[2])

        r_kalman_filter.priori()
        r_kalman_filter.update(*gyro, *accel)

        current_state = kalman_filter.get_state()
        current_state_r = r_kalman_filter.get_state()

        apogee_est = apogee_estimator.predict_apogee(current_state[0:3])

        flap_ext = controller.get_flap_extension(time_stamp > prop.delay and np.linalg.norm(motor.get_thrust(time_stamp)) <= 0, apogee_est)

        # flap_ext will be passed by kalman filter
        rocket.set_motor_mass(time_stamp)
        # rocket.motor_mass = motor.get_mass(time_stamp)

        x, alpha = sim.RK4(x, dt, time_stamp, flap_ext)
        time_stamp += dt

        addToDict(x, baro_alt, accel, bno_ang_pos, gyro, current_state, current_state_r, alpha, apogee_est, rocket.rocket_total_mass, rocket.motor_mass, flap_ext)

    t_end = time.time() - t_start
    print(f"Time: {t_end:.2f}")


if __name__ == '__main__':
    x0 = np.zeros((6, 3))
    x0[3] = [0, 0.05, 0]
    dt = 0.01
    simulator(x0, dt)

    # plot entries in sim_dict
    # print(np.array(sim_dict["pos"]))
    # # print(sim_dict["time"])
    # plt.plot(sim_dict["time"], np.array(sim_dict["pos"])[:,0])
    # plt.plot(sim_dict["time"], np.array(sim_dict["vel"])[:,0])
    # plt.show()
    # plotter.plotter(sim_dict=sim_dict)

    print("Writing to file...")

    record = []
    for point in range(len(sim_dict["time"])):
        cur_point = []
        cur_point.append(str(sim_dict["time"][point]))
        cur_point += list(map(str, sim_dict["pos"][point]))
        cur_point += list(map(str, sim_dict["vel"][point]))
        cur_point += list(map(str, sim_dict["accel"][point]))
        cur_point += list(map(str, sim_dict["ang_pos"][point]))
        cur_point += list(map(str, sim_dict["ang_vel"][point]))
        cur_point += list(map(str, sim_dict["ang_accel"][point]))
        cur_point += map(str, list([sim_dict["alpha"][point]]))
        cur_point += map(str, list([sim_dict["rocket_total_mass"][point]]))
        cur_point += map(str, list([sim_dict["motor_mass"][point]]))
        cur_point += map(str, list([sim_dict["flap_ext"][point]]))
        cur_point += map(str, list([sensor_dict["baro_alt"][point]]))
        cur_point += map(str, list([sensor_dict["imu_accel_x"][point]]))
        cur_point += map(str, list([sensor_dict["imu_accel_y"][point]]))
        cur_point += map(str, list([sensor_dict["imu_accel_z"][point]]))
        cur_point += map(str, list([sensor_dict["imu_ang_pos_x"][point]]))
        cur_point += map(str, list([sensor_dict["imu_ang_pos_y"][point]]))
        cur_point += map(str, list([sensor_dict["imu_ang_pos_z"][point]]))
        cur_point += map(str, list([sensor_dict["imu_gyro_x"][point]]))
        cur_point += map(str, list([sensor_dict["imu_gyro_y"][point]]))
        cur_point += map(str, list([sensor_dict["imu_gyro_z"][point]]))
        cur_point += map(str, list([sensor_dict["apogee_estimate"][point]]))
        cur_point += map(str, list(kalman_dict["x"][point]))
        cur_point += map(str, list(kalman_dict["y"][point]))
        cur_point += map(str, list(kalman_dict["z"][point]))
        cur_point += map(str, list(kalman_dict["rx"][point]))
        cur_point += map(str, list(kalman_dict["ry"][point]))
        cur_point += map(str, list(kalman_dict["rz"][point]))

        record.append(cur_point)

    output_file = os.path.join(os.path.dirname(__file__), prop.output_file)
    with open(output_file, 'w') as f:
        f.write("time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,accel_x,accel_y,accel_z,ang_pos_x,ang_pos_y,ang_pos_z,ang_vel_x,ang_vel_y,ang_vel_z,ang_accel_x,ang_accel_y,ang_accel_z,alpha,rocket_total_mass,motor_mass,flap_ext,baro_alt,imu_accel_x,imu_accel_y,imu_accel_z,imu_ang_pos_x,imu_ang_pos_y,imu_ang_pos_z,imu_gyro_x,imu_gyro_y,imu_gyro_z,apogee_estimate,kalman_pos_x,kalman_vel_x,kalman_accel_x,kalman_pos_y,kalman_vel_y,kalman_accel_y,kalman_pos_z,kalman_vel_z,kalman_accel_z,kalman_rpos_x,kalman_rvel_x,kalman_raccel_x,kalman_rpos_y,kalman_rvel_y,kalman_raccel_y,kalman_rpos_z,kalman_rvel_z,kalman_raccel_z\n")
        for point in record:
            f.write(f"{','.join(point)}\n")
