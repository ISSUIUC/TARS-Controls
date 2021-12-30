import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix
from matplotlib import pyplot as plt
from statistics import mean

class density:
    def density_func(z):
        #? Defining a few constants 
        #* temperature under standard condition (15 degrees C at sealevel) kelvin
        T_0 = 288.16 
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
        rho = rho_0*(1 - ((b*z)/T_0))**(g/(R*b))*(T_0/(T_0 - b*z))
        
        return rho

class conversion:
    def ft_to_m(measurement): 
        return (measurement / 3.2808) 

    def c_to_k(measurement):
        return (measurement + 273.15)

    def deg_to_rad(measurement):
        return (measurement*np.pi/180)

    def rad_to_deg(measurement):
        return (measurement*180/np.pi)

class speed_sound:
    def speed_sound(altitude):
        #* the first 11000 m of flight, the temperature varies linearly with altitude (k/m)
        dt_dh = 6.5e-3 
        #* specific heat ratio of air 
        r = 1.4 
        #* ideal gas constant of air R 
        R = 287.05
        #* temperature under standard condition (15 degrees C at sealevel) kelvin
        T_0 = 288.16 
        #* first derivative 
        da_dh = 0.5*dt_dh*np.sqrt(r*R/(T_0 + dt_dh*altitude))
        #* linearization 
        a_H = a_h + da_dh*(altitude - h)

        return a_H