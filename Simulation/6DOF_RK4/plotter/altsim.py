# File for onbaording sim
import numpy as np
import time


class Sim:
    def __init__(self, mass, thrust, cDrag, crossA, dm, dt, x0, temp0, pressure0):
        self.mass = mass
        self.thrust = thrust
        self.cDrag = cDrag
        self.crossA = crossA
        self.dm = dm
        self.dt = dt
        self.temp0 = temp0
        self.pressure0 = pressure0
        self.time = 0
        self.x = x0
        self.v = 0
        self.G = 6.6743 / np.power(10,11)
        self.mEarth = 5.97219 * np.power(10,24)
        self.rEarth = 6.378137 * np.power(10,6)
        self.Fg = self.mass * self.G * self.mEarth / np.power(self.rEarth+self.x,2)
        self.fNet = 0

    def position(self):
        self.v = self.fNet*self.dt
        self.x += self.v*self.dt
        
    def aDrag(self):
        self.temp = self.temp0 - .0065*self.x
        self.pressure = self.pressure0*((self.temp0/self.temp)**((9.81*.029)/(8.314*.0065)))
        self.density = (.029*self.temp)/(8.314*self.pressure)
        self.drag = .5*self.density*self.v*self.v*self.cDrag*self.crossA

    def forceGravity(self):
        self.mass += self.dm * self.dt
        self.g = self.G * self.mEarth / (self.rEarth+self.x ** 2)
        self.Fg = self.mass * self.g

    def netForce(self):
        self.fNet = self.thrust - self.drag - self.Fg

    def update(self):
        self.position()
        self.aDrag()
        self.forceGravity()
        self.netForce()

if __name__ == '__main__':
    start = np.array([3,2,1,"GO"])
    for i in range(0,4):
        print(start[i])
        time.sleep(1)

    launch = Sim(1, 6700, .75, 4, 5, .01, 0, 33, 102500)
    while launch.position >= 0:
        launch.update()
    