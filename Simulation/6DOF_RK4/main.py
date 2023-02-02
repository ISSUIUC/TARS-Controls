import numpy as np
import matplotlib.pyplot as plt
import os

import motor
import properties as prop
import simulator as sim
import plotSIM as plotter
import sensors
import time

motor = motor.Motor()

sim_dict = {
    "pos":[],
    "vel": [],
    "accel": [],
    "ang_pos":[],
    "ang_vel": [],
    "ang_accel": [],
    "alpha": [],
    "time": []
    }

sensor_dict ={
    "baro_alt":[],
    "accel_x": [],
    "accel_y": [],
    "accel_z":[],
    "bno_ang_pos_x":[],
    "bno_ang_pos_y":[],
    "bno_ang_pos_z":[],
    "gyro_x":[],
    "gyro_y": [],
    "gyro_z": []
}

def simulator(x0, dt) -> None:
    '''
    Method which handles running the simulation and logging sim data to dict

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
    time_stamp = 0
    idle_time = 0 # time in seconds before launch
    while time_stamp < idle_time:
        time_stamp += dt
    
    print("Ignition")

    motor.ignite(time_stamp*dt)
    # # while x[1][prop.vertical] > prop.apogee_thresh and x[0][prop.vertical] > prop.start_thresh:
    start = True
    t_start = time.time()
    while x[1,0] >= 0 or start:
        if start:
            start = False
            
        # Get sensor data
        baro_alt = sensors.get_barometer_data(x)
        accel_x, accel_y ,accel_z = sensors.get_accelerometer_data(x)
        gyro_x, gyro_y, gyro_z = sensors.get_gyro_data(x)
        bno_ang_pos_x, bno_ang_pos_y, bno_ang_pos_z = sensors.get_bno_orientation(x)
        
        # Append to sensor_dict
        sensor_dict["baro_alt"].append(baro_alt)
        sensor_dict["accel_x"].append(accel_x)
        sensor_dict["accel_y"].append(accel_y)
        sensor_dict["accel_z"].append(accel_z)
        sensor_dict["bno_ang_pos_x"].append(bno_ang_pos_x)
        sensor_dict["bno_ang_pos_y"].append(bno_ang_pos_y)
        sensor_dict["bno_ang_pos_z"].append(bno_ang_pos_z)
        sensor_dict["gyro_x"].append(gyro_x)
        sensor_dict["gyro_y"].append(gyro_y)
        sensor_dict["gyro_z"].append(gyro_z)
        # Kalman Filter stuff goes here
        
        # flap_ext will be passed by kalman filter
        prop.motor_mass = motor.get_mass(time_stamp)

        # Update Simulator Log
        sim_dict["pos"].append(x[0])
        sim_dict["vel"].append(x[1])
        sim_dict["accel"].append(x[2])
        sim_dict["ang_pos"].append(x[3])
        sim_dict["ang_vel"].append(x[4])
        sim_dict["ang_accel"].append(x[5])
        sim_dict["time"].append(sim_dict["time"][-1]+dt if len(sim_dict["time"]) > 0 else 0)

        x,alpha = sim.RK4(x, dt, time_stamp)
        sim_dict["alpha"].append(alpha)
        time_stamp += dt
    t_end = time.time() - t_start
    print("Time: ", t_end)


if __name__ == '__main__':
    x0 = np.zeros((6,3))
    x0[3] = [0, .05, 0]
    x0[1] = [0,0,0]
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
        record.append(cur_point)
    
    output_file = os.path.join(os.path.dirname(__file__), prop.output_file)
    with open(output_file, 'w') as f:
        f.write("time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,accel_x,accel_y,accel_z,ang_pos_x,ang_pos_y,ang_pos_z,ang_vel_x,ang_vel_y,ang_vel_z,ang_accel_x,ang_accel_y,ang_accel_z,alpha\n")
        for point in record:
            f.write(f"{','.join(point)}\n")