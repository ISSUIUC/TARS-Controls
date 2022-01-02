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
    def speed_sound(altitude,a_h,h):
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
        minv = np.min(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
        m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
        m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
        min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
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
        minv = np.min(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
        m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
        m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
        min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
        ax.legend(handles=[max,m1,m2,m3,min])
        return

class inertia:
    def I_new(m,r_m):
        #m is added mass (in grams) to nosecone 
        # r_m is the distance (in cm) from the tip of the nosecone to the CG of added mass (can be estimated) 
    
        m = input("dry mass addition(kg):");m = float(m)
        r_m = input("dry mass location(m):");r_m = float(r_m)

    
        #constants obtained from open rocket data
        r_CG_0 = 167.67 #cm #at burnout 
        m0 = 21221 #g #at burnout
        Ixx_0 = 0.030245 #kg*m^2
        Iyy_0 = 15.841 #kg*m^2
        Izz_0 = 15.841 #kg*m^2
        
        #unit conversions
        r_CG_0 = r_CG_0/100 #converted to m
        m0 = m0/1000 #converted to kg
        
        #finding new center of mass location
        new_m = m0 + m
        new_r_CG = (m0*r_CG_0 + m*r_m) / new_m
        d1 = new_r_CG - r_m #distance between new CG and the location where the mass is added
        d2 = r_CG_0 - new_r_CG #distance between new and old CGs
        
        #calculating new moments of inertia around new CG 
        Ixx_new = Ixx_0
        Iyy_new = Iyy_0 + m*d1**2 + m0*d2**2
        Izz_new = Iyy_new
        return(Ixx_new,Iyy_new,Izz_new)

class rotation: 
    
    def yaw(phe):
        #* rotation about the z-axis, takes in angle phe in the function (rad)
        R_y = np.array([[np.cos(phe), -np.sin(phe), 0],[np.sin(phe), np.cos(phe), 0],[0, 0, 1]])
        return R_y

    def pitch(theta):
        #* rotation about the y-axis, takes in angle theta in the function (rad)
        R_p = np.array([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0],[-np.sin(theta), 0, np.cos(theta)]])
        return R_p
    
    def roll(phi):
        #* rotation about the x-axis, take sin angle phi in the function (rad)
        R_r = np.array([[1, 0, 0],[0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
        return R_r

    def body_aero(Vx_b, Vy_b, Vz_b):
        #* rotation matrix that transforms the body frame to the aerodyanmic frame 
        #* takes in beta (side-slip angle), and alpha (angles of attack)
        #* needs the velocities in the body frame to calculate beta and alpha
        alpha = np.arctan(Vz_b/Vx_b)
        beta = np.arctan(Vy_b/(np.sqrt(Vx_b**2 + Vz_b**2)))
        R_ba = np.array([[np.cos(beta)*np.cos(alpha), np.sin(beta), np.cos(beta)*np.sin(alpha)], [-np.sin(beta)*np.cos(alpha), np.cos(beta), 
        -np.sin(beta)*np.sin(alpha)],[-np.sin(alpha), 0, np.cos(alpha)]])
        return R_ba
    