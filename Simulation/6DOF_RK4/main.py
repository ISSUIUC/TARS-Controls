import numpy as np
import matplotlib.pyplot as plt

import motor
import properties as prop
import simulator as sim
import plotSIM as plotter

motor = motor.Motor()

sim_dict = {
    "pos":[],
    "vel": [],
    "accel": [],
    "time": []
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
    while x[1,0] > 0 or start:
        if start:
            start = False
        # Kalman Filter stuff goes here
        # flap_ext will be passed by kalman filter
        prop.motor_mass = motor.get_mass(time_stamp)

        x = sim.newtonProp(x, dt, time_stamp)

        # Update Simulator Log
        sim_dict["pos"].append(x[0])
        sim_dict["vel"].append(x[1])
        sim_dict["accel"].append(x[2])
        sim_dict["time"].append(sim_dict["time"][-1]+dt if len(sim_dict["time"]) > 0 else 0)

        time_stamp += dt


if __name__ == '__main__':
    x0 = np.zeros((6,3))
    dt = 0.05
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
        record.append(cur_point)
    
    with open(prop.output_file, 'w') as f:
        f.write("time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,accel_x,accel_y,accel_z\n")
        for point in record:
            f.write(f"{','.join(point)}\n")