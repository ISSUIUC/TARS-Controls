import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix
from matplotlib import pyplot as plt
from statistics import mean


def density_func(z):
    #? Defining a few constants 
    #* temperature under standard condition (15 degrees C at sealevel) kelvin
    T_o = 288.16 
    #* pressure under standard condition in (Pa)
    P_0 = 101325
    #* Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
    b = 0.0065
    #* gravitational constant 
    g = 9.81
    #* air density under standard condition (kg/m3)
    rho_0 = 1.225 
    #* ideal gas constant (J/kg K)
    R = 287.05

    #? returning the density at an altitude Z
    rho = rho_0*(1 - ((b*z)/T_o))**(g/(R*b))*(T_o/(T_o - b*z))
    
    return rho

h = np.arange(0,10000,0.1)
RHO = []
for num in h:
    RHO.append(density_func(num))

fig = plt.figure(dpi = 200)
plt.plot(h, RHO, label="Air Density (kg/m3)"); plt.ylabel("Density");plt.xlabel("Altitude (m)")
plt.show()