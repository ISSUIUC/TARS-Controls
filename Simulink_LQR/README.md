# IREC LQR Simulation and Calculation

The goal is to use LQR to optimize our actuator inputs to the 4 servos, which change the exposed surface area of the flaps to control roll and drag.
The 4 servos come in two pairs. The first pair are diametrically opposite each other and are angled angled at 45°. The other pair is angled at -45° and is offset by 90 degrees from the previous pair, so the four servos are located at 4 equally spaced points.

## Equations of Motion
The following equations of motion were used to describe the dynamics of the system.

![EOM](SimulinkLQR/images/EOM.jpg)

## Variables
The descriptions of the constants, states, inputs and the outputs of the system are shown below. The LQR is setup and the gains are calculated according to the below ordering of variables.

![vars](.../images/Vars.jpg)

## Cost function
The cost function below is used in the implementation of the LQR system.

![cost](.../images/costf.jpg)

## Simulation
The silumation is done in the Simulink model "LQRsim.slx". The equations of motion above were modeled in a Simulink block called 'Plant'. The input vector u enters the plant and the output y is the result. The output is reordered to make x, multiplied by -K (the LQR gains) to get the input u which is fed back into the Plant. A saturation limit and time delay are also included in the model for realistic results. The scopes are for seeing the results.

![sim](.../images/sim.jpg)

## Calculation
All the variables and calculation of the LQR gains is done in steps in the LQRcalc.m script.
