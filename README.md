# TARS-CONTROLS
Illinois Space Society's Large Rocketry Guidance, Navigation, and Controls Team simulation tools repository

We test various control design methodologies here with python before implementing them in SILSIM and flight software.

In order to run these programs, we recommend installing **Anaconda**. Run the command: `conda env create -f ISS-GNC.yaml` in your command line to create a custom conda environment and install the required packages.

The `Legacy/Control_Design` folder contains the python implementation of the PID roll controller that was used on the **Endurance** launch vehicle in October 2021.

`Simulation/Legacy/1dof_sim_RK4.py` is an RK4 based 1-DOF python simulation that was used to develop the active drag control system and extended kalman filter that is implemented on the **Intrepid** class of launch vehicles in 2021-2022. 

Our current simulation: `Simulation/6DOF_RK4` which can be run from `6DOF_RK4/simulation.main` was used to develop the Kalman Filter and controller that was used on **Intrepid-III**, which placed **2nd in the 30K COTS** category at **Spaceport America Cup 2023**. This simulation is being updated to include a model for rocket descent, monte-carlo configuration, and support for staged rockets. This will be used to develop flight event logic such as stage separation, air-lit motor ignition, parachute deployment and parachute disreefing.


