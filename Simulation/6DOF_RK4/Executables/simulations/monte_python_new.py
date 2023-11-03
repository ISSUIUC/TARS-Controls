import numpy as np
import os
import sys
import shutil
import time

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..','..','Core')))

import simulation as simulation
import dynamics.rocket as rocket_model
import dynamics.propagator as propagtor
import environment.atmosphere as atmosphere
import util.data_process as data_process

import properties as prop
import data_loader as dataloader
import pysim as pysim

def monte(x0, sim, rocket, atm, dt, kf_dt, samples, output_folder, target_size, clear_contents=True) -> None:
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
    if clear_contents and os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        os.mkdir(f"{output_folder}/figures")
        os.mkdir(f"{output_folder}/SimData")

    print("Running Monte Python")
    t_start = time.time()
    # Nominal run
    print("Calculating nominal trajectory")
    record = pysim.pysim(dt, kf_dt)
    data_process.save_record(record, os.path.join(output_folder, "SimData", "nominal.csv"), target_size=target_size)

    # Monte Carlo runs
    for i in range(samples):
        print(f"Running sample {i+1} of {samples}")
        pysim.simulator(x0, sim, rocket, atm, dt, kf_dt)
        record = rocket.to_list()
        data_process.save_record(record, os.path.join(output_folder, "SimData", f"sim_data_{i}.csv"), target_size=target_size)
    t_end = time.time() - t_start
    print(f"Monte Python Runtime: {t_end:.2f} seconds")


config = dataloader.config
x0 = np.zeros((6,3))
x0[3] = np.array([0.0, 0.05, 0.01])
dt = 0.1
kf_dt = 0.1

target_size = 880
samples = 2

atm = atmosphere.Atmosphere(wind_direction_variance_stddev=0.2, 
                            wind_magnitude_variance_stddev=3.0,
                            enable_direction_variance=True,
                            enable_magnitude_variance=True, 
                            nominal_wind_magnitude=2,
                            nominal_wind_direction=np.array([0,1,0]))
stages = []
for stage in config['rocket']['stages'][1:]:
    stages.append(rocket_model.Rocket(stage, atm=atm))
rocket = rocket_model.Rocket(config['rocket']['stages'][0], atm=atm, stages=stages)

sim = propagtor.Simulator(atm=atm, rocket=rocket)
output_folder = os.path.join(os.path.dirname(__file__), config["meta"]["monte_carlo_output_folder"])

monte(x0, sim, rocket, atm, dt, kf_dt, samples, output_folder, target_size)