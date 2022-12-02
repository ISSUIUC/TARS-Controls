import properties as prop
import simulator as sim
sim_dict = {
    "pos":[],
    "vel": [],
    "accel": [],
    "time": [0]
    }

def simulator(x0, dt):
    '''
    x0 --> initial state
    dt --> time step
    '''
    x = x0.copy()
    while x[1][prop.vertical] > prop.apogee_thresh and x[0][prop.vertical] > prop.start_thresh:
        # Kalman Filter stuff goes here
        # flap_ext will be passed by kalman filter

        x = sim.RK4(x, dt)

        # Update Simulator Log
        sim_dict["pos"].append(x[0])
        sim_dict["vel"].append(x[1])
        sim_dict["accel"].append(x[2])
        sim_dict["time"].append(sim_dict["time"][-1]+dt)