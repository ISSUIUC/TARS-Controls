import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import dynamics.motor as motor
import dynamics.forces as forces

class Rocket:

    # Center of mass of entire body
    cm = np.array([3.34-2.31, 0., 0.])
    cp = np.array([3.34-2.71, 0., 0.])
    # Center of mass of rocket without motor
    cm_rocket = np.array([3.34-1.86, 0., 0.])
    # Center of mass of the motor
    cm_motor = np.array([0.3755, 0., 0.])

    # Motor Properties M2500:
    impulse = 9671.0  # Ns
    motor_mass = 8.064  # Kg
    delay = 60  # s
    motor_lookup_file = '../lookup/m2500.csv'

    # # Motor Properties N5800:
    # impulse = 20145.7 # Ns
    # motor_mass = 14.826 # Kg
    # delay = 0 # s
    # motor_lookup_file = '../lookup/n5800.csv'

    # # Motor Properties N2540:
    # impulse = 17907 # Ns
    # motor_mass = 10.700 # Kg
    # delay = 0 # s
    # motor_lookup_file = '../lookup/n2540.csv'

    # rocket mass w/out motor
    rocket_dry_mass = 14.691
    # rocket mass with motor
    rocket_total_mass = rocket_dry_mass + motor_mass
    # radius of rocket
    r_r = 0.0508
    # length of rocket
    l = 3.34
    # area w/out flaps
    A = math.pi*r_r**2
    # side profile area
    A_s = 2*r_r*l
    # flap max estension length (m)
    max_ext_length = .0178

    motor = None
    forces = None
    
    def __init__(self, cm_rocket=np.array([3.34-1.86, 0., 0.]), cm_motor=np.array([0.3755, 0., 0.]),impulse=9671.0,
                 motor_mass=8.064, delay=60, motor_lookup_file='../lookup/m2500.csv',rocket_dry_mass=14.691,r_r=0.0508,l=3.34,
                 max_ext_length=0.0178, atm=None):
        self.cm_rocket = cm_rocket
        self.cm_motor = cm_motor

        self.impulse = impulse
        self.motor_mass = motor_mass
        self.delay = delay
        self.motor_lookup_file = motor_lookup_file

        self.rocket_dry_mass = rocket_dry_mass
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        self.r_r = r_r
        self.l = l
        self.A = math.pi * self.r_r ** 2
        self.A_s = 2 * self.r_r * self.l
        self.max_ext_length = max_ext_length
        self.atm = atm
        
        # need to change motor and forces constructors to refer to properties from this class rather than properties file
        self.motor = motor.Motor(self.rocket_dry_mass, 
                                 self.cm, 
                                 self.cm_rocket, 
                                 self.cm_motor, 
                                 self.rocket_dry_mass, 
                                 impulse=self.impulse, 
                                 mass=self.motor_mass, 
                                 delay=self.delay, 
                                 lookup_file=self.motor_lookup_file)
        self.forces = forces.Forces(self.max_ext_length,
                                    self.cm,
                                    self.cp,
                                    self.A,
                                    self.A_s,
                                    self.rocket_dry_mass,
                                    self.motor,
                                    self.atm)

    def set_motor_mass(self, timestamp):
        """Sets the mass of the motor at a given time
        
        Args:
            timestamp (float): Time in seconds
        """
        self.motor_mass = self.motor.get_mass(timestamp)
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        
    def I(self, total_mass=rocket_total_mass): 
        """Returns the inertia matrix of the rocket
        
        Args:
            total_mass (float, optional): Total mass of the rocket. Defaults to rocket_total_mass.
        
        Returns:
            np.array: Inertia matrix of the rocket
        """
        return np.diag([(1/2) * total_mass * self.r_r**2,
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2),
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2)])

    
    def I_inv(self, total_mass=rocket_total_mass): 
        """Returns the inverse of the inertia matrix of the rocket
        
        Args:
            total_mass (float, optional): Total mass of the rocket. Defaults to rocket_total_mass.
        """
        return np.diag([1/((1/2) * total_mass * self.r_r**2),
                                       1/((total_mass/12) * (self.l**2 + 3*self.r_r**2)),
                                       1/((total_mass/12) * (self.l**2 + 3*self.r_r**2))])
    

# if __name__ == '__main__':
    # rocket = Rocket(motor_mass=5,rocket_dry_mass=2)
    # print(rocket.rocket_total_mass)
