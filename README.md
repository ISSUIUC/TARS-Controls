# TARS-CONTROLS
Illinois Space Society's Large Rocketry Guidance, Navigation, and Controls Team simulation tools and tutorials repository

We test various control design methodologies here with python before implementing them in SILSIM and flight software.

In order to run these programs, we recommend installing **Anaconda** and using the `-conda install <package_name>` to install the required packages.

The `Control_Design` folder contains the python implementation of the PID roll controller that was used on the **Endurance** launch vehicle in October 2021.

`1dof_sim_RK4.py` is an RK4 based 1-DOF python simulation that was used to develop the active drag control system and extended kalman filter that is implemented on the **Intrepid** class of launch vehicles. The Kalman filter that is being developed here and eventually in `6dof_sim.py` will also be implemented for testing on individual's rockets that are in the ISS **L1 Research Fleet** for test and evaluation of new components and algorithms. `1dof_sim_RK4.py` will be used as a framework to improve `6dof_sim.py`.


