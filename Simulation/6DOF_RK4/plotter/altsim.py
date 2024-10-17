# File for onbaording sim
import matplotlib.pyplot as plt
import numpy as np
import time


class Sim:
    def __init__(self, mass, thrust, cDrag, crossA, motorMass, dm, dt, x0, temp0, pressure0):
        self.mass = mass
        self.thrust = thrust
        self.cDrag = cDrag
        self.crossA = crossA
        self.dm = dm
        self.dt = dt
        self.temp0 = temp0 + 274.15
        self.pressure0 = pressure0
        self.motorMass = motorMass
        self.time = 0
        self.x = x0
        self.v = 0
        self.a = 0
        self.G = 6.6743/(10**11)
        self.mEarth = 5.97219*(10**24)
        self.rEarth = 6.378137*(10**6)
        self.fNet = 0

    def position(self):
        self.a = self.fNet/self.mass
        self.v += self.a*self.dt
        self.x += self.v*self.dt
        
    def aDrag(self):
        self.temp = self.temp0 - .0065*self.x
        if(self.temp > 0):
            self.pressure = self.pressure0*((self.temp0/self.temp)**((9.81*.029)/(8.314*.0065)))
        else:
            self.pressure = -self.pressure0*((self.temp0/(-self.temp))**((9.81*.029)/(8.314*.0065)))
        self.density = (.029*self.temp)/(8.314*self.pressure)
        self.drag = .5*self.density*self.v*self.v*self.cDrag*self.crossA

    def forceGravity(self):
        self.g = (self.G * self.mEarth)/((self.rEarth+self.x)**2)
        self.Fg = self.mass * self.g

    def endStage(self):
        if(self.motorMass <= 0):
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
        self.endStage()
        self.netForce()
        self.position()
        self.time += self.dt

if __name__ == '__main__':
    launch = Sim(50, 2000, .75, .25, 40, 1, .01, 0, 33, 102500)
    x = [0]
    t = [0]
    drag = [0]
    while launch.x >= 0:
        launch.update()
        x.append(launch.x)
        t.append(launch.time)
        drag.append(launch.drag)
    plt.plot(t, x)
    plt.xlabel('time(s)')
    plt.ylabel('height(m)')
    plt.show()
    # plt.plot(t, drag)
    # plt.xlabel('time(s)')
    # plt.ylabel('drag(N)')
    # plt.show()
    