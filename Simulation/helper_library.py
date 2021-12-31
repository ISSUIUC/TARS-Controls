import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix
from matplotlib import pyplot as plt
from statistics import mean
from mpl_toolkits import mplot3d
import matplotlib.pyplot as pyplt
import matplotlib.patches as mpatches
import numpy as np 

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


class plot:
    # plots a trajectory based on an array of 3d points
    def plot_3d(points):
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = points[:, 0]
        x = points[:, 1]
        y = points[:, 2]
        ax.scatter(x,y,z, c=[.75,0,0], alpha=1)
        return

    # plots a trajectory based on an array of 3d points 
    # color maps the velocity based on a vector of velocity magnitudes
    def plot_3d_vel(points, velocities):
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = points[:, 0]
        x = points[:, 1]
        y = points[:, 2]
        maxv = np.max(velocities)
        color = np.zeros((velocities.shape[0], 3))
        med = np.median(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[.75,0,0])
        m2 = mpatches.Patch(color=[.5,0,0], label=med)
        m3 = mpatches.Patch(color=[.25,0,0])
        min = mpatches.Patch(color=[0,0,0], label= np.min(velocities))
        ax.legend(handles=[max,m1,m2,m3,min])
        return

    # plots trajectory with velocity estimated based on a timestep
    def plot_3d_est(points, timestep):
        velocities = np.zeros(points.shape[0])
        velocities[0] = abs(np.linalg.norm(points[1] - points[0]))/timestep
        for i in range(1, velocities.shape[0]):
            velocities[i] = abs(np.linalg.norm(points[i] - points[i-1]))/timestep
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = points[:, 0]
        x = points[:, 1]
        y = points[:, 2]
        maxv = np.max(velocities)
        color = np.zeros((velocities.shape[0], 3))
        med = np.median(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[.75,0,0])
        m2 = mpatches.Patch(color=[.5,0,0], label=med)
        m3 = mpatches.Patch(color=[.25,0,0])
        min = mpatches.Patch(color=[0,0,0], label= np.min(velocities))
        ax.legend(handles=[max,m1,m2,m3,min])
        return
