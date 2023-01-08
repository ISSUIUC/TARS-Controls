import numpy as np
import matplotlib.pyplot as plt

import motor
import properties as prop
import simulator as sim

motor = motor.Motor()

sim_dict = {
    "pos":[],
    "vel": [],
    "accel": [],
    "time": []
    }

def simulator(x0, dt):
    '''
    x0 --> initial state
    dt --> time step
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
    while x[0][1] > 0 or start:
        if start:
            start = False
        # Kalman Filter stuff goes here
        # flap_ext will be passed by kalman filter
        prop.motor_mass = motor.get_mass(time_stamp)

        x = sim.RK4(x, dt, time_stamp)

        # Update Simulator Log
        sim_dict["pos"].append(x[0])
        sim_dict["vel"].append(x[1])
        sim_dict["accel"].append(x[2])
        sim_dict["time"].append(sim_dict["time"][-1]+dt if len(sim_dict["time"]) > 0 else 0)

        time_stamp += dt


if __name__ == '__main__':
    x0 = np.zeros((3,6))
    dt = 0.01
    simulator(x0, dt)
    # plot entries in sim_dict
    print(np.array(sim_dict["pos"])[:,0])
    # print(sim_dict["time"])
    plt.plot(sim_dict["time"], np.array(sim_dict["pos"])[:,0])
    plt.show()
