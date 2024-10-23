#Link to equations used: <a href="https://pages.vassar.edu/magnes/2019/05/12/computational-simulation-of-rocket-trajectories/"> </a>
# File for onbaording sim
import matplotlib.pyplot as plt
import numpy as np
import time


class Sim:
    def __init__(self, mass, thrust, cDrag, radius, motorMass, dm, dt, x0, temp0, pressure0):
        self.mass = mass # Mass of the rocket in kg
        self.thrust = thrust # Thrust of the rocket in N
        self.cDrag = cDrag # Coefficient of drag
        self.crossA = np.pi*radius*radius # Cross sectional area of the rocket in m from radius of the rocket (pi*(r^2)
        self.dm = dm # Change in mass with respect to time as motor fires in kg/s
        self.dt = dt # How much you want to step time for each calculation
        self.temp0 = temp0 + 274.15 # ground temperature in C (is converted to K in this line)
        self.pressure0 = pressure0 # ground pressure in pascals
        self.motorMass = motorMass # mass of the motor
        self.time = 0 # Time (s) after launch
        self.x = x0 # Height (m) of the rocket above sea level
        self.v = 0 # Velocity (m/s) of the rocket
        self.a = 0 # Acceleration (m/s^2) of the rocket
        self.fNet = 0 # Net force (N) acting on the rocket: combination of thrust, drag and gravity
        self.g = 9.81 # Acceleration (m/s^2) due to gravity

    def position(self):
        self.a = self.fNet/self.mass # F = ma -> a = F/m
        self.v += self.a*self.dt # v = a*dt
        self.x += self.v*self.dt # x = v*dt
        
    def aDrag(self):
        L = .0065 # L is the temperature lapse rate (estimated as .0065)
        self.temp = self.temp0 - self.L*self.x # *Temperature(x)* T(x) = temp0 - Lx
        M = .029 # M is the molar mass of air 0.029 kg/mol
        R = 8.314 # R is the gas constant 8.314 J*(K^-1)*(mol^-1)
        if(self.temp > 0):
            self.pressure = self.pressure0*((self.temp0/self.temp)**((self.g*M)/(R*L)))
            # *Pressure(x)* P(x) = pressure0*(temp0/T(x))^((g*M)/(RL))
        else:
            self.pressure = -self.pressure0*((self.temp0/(-self.temp))**((self.g*M)/(R*L)))
            # same as above but prevents imaginaries when T(x) is negative
        self.density = (M*self.temp)/(R*self.pressure) # *AirDensity(x)* p(x) = (M*T(x))/(R*P(x))
        self.drag = .5*self.density*self.v*self.v*self.cDrag*self.crossA # Drag = .5*p(x)*(v^2)*DragCoefficient*CrossSectionalArea

    def forceGravity(self):
        self.Fg = self.mass * self.g # Fg = mg

    def updateMass(self): 
        if(self.motorMass <= 0): # out of fuel
            self.thrust = 0
            self.dm = 0
        self.mass -= self.dm*self.dt 
        self.motorMass -= self.dm*self.dt

    def netForce(self):
        if(self.v >= 0):
            self.fNet = self.thrust - self.drag - self.Fg
        else:
            self.fNet = self.drag - self.Fg

    def update(self):
        self.aDrag()
        self.forceGravity()
        self.updateMass()
        self.netForce()
        self.position()
        self.time += self.dt

if __name__ == '__main__':
    launch = Sim(50, 2000, .75, .25, 40, 1, .01, 0, 33, 102500) # Values we inputted are made up
    x = [0] # array to store altitudes to plot
    t = [0] # array to store time to plot
    while launch.x >= 0:
        launch.update()
        x.append(launch.x)
        t.append(launch.time)
    plt.plot(t, x)
    plt.xlabel('time(s)')
    plt.ylabel('height(m)')
    plt.show()
    