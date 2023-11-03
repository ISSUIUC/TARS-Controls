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
import os
import sys
import time

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..','..','Core')))

import simulation as simulation
import dynamics.rocket as rocket_model
import dynamics.propagator as propagtor
import environment.atmosphere as atmosphere
import util.data_process as data_process

import data_loader as dataloader

def simulator(x0, sim, rocket, atm, dt, kf_dt) -> None:
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
    simulator = simulation.Simulation(sim, rocket, atm, dt, x0, kf_dt, stages=rocket.stages)
    simulator.idle_stage()
    t_start = time.time()
    simulator.run_stages()
    simulator.coast()
    t_end = time.time() - t_start
    print(f"Runtime: {t_end:.2f} seconds")

def pysim(dt=0.01, kf_dt=0.01, config=dataloader.config):
    # Load desired config file
    x0 = np.zeros((6, 3))
    x0[3] = [0, 0.05, 0]

    atm = atmosphere.Atmosphere()

    stages = []
    for stage in config['rocket']['stages'][1:]:
        stages.append(rocket_model.Rocket(stage, atm=atm))
    rocket = rocket_model.Rocket(config['rocket']['stages'][0], atm=atm, stages=stages)

    sim = propagtor.Simulator(atm=atm, rocket=rocket)

    simulator(x0, sim, rocket, atm, dt, kf_dt)
    return rocket.to_list()

if __name__ == "__main__":
    config = dataloader.config
    record = pysim()
    print("Simulation complete")
    output_file = os.path.join(os.path.dirname(__file__), config["meta"]["output_file"])
    data_process.save_record(record, output_file)