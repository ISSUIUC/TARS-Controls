import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import dynamics.motor as motor
import dynamics.forces as forces

class Rocket:
    motor = None
    forces = None
    stage_config = None
    
    def __init__(self, stage_config, atm=None):
        self.stage_config = stage_config
        self.cm_rocket = stage_config["rocket_body"]["structure_cm"]
        self.cm_motor = stage_config["motor"]["cm"]
        self.cm = stage_config["rocket_body"]["combined_cm"]
        self.cp = stage_config["rocket_body"]["combined_cp"]

        self.impulse = stage_config["motor"]["impulse"]
        self.motor_mass = stage_config["motor"]["motor_mass"]
        self.delay = stage_config["motor"]["delay"]
        self.motor_lookup_file = stage_config["motor"]["motor_lookup_file"]

        self.rocket_dry_mass = stage_config["rocket_body"]["dry_mass"]
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        self.r_r = stage_config["rocket_body"]["radius"]
        self.l = stage_config["rocket_body"]["length"]
        self.A = math.pi * self.r_r ** 2
        self.A_s = 2 * self.r_r * self.l
        self.max_ext_length = stage_config["flaps"]["max_ext_length"]
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
                                    stage_config["rocket_body"]["rasaero_lookup_file"],
                                    self.atm)

    def set_motor_mass(self, timestamp):
        """Sets the mass of the motor at a given time
        
        Args:
            timestamp (float): Time in seconds
        """
        self.motor_mass = self.motor.get_mass(timestamp)
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
        
    def I(self, total_mass): 
        """Returns the inertia matrix of the rocket
        
        Args:
            total_mass (float, optional): Total mass of the rocket. Defaults to rocket_total_mass.
        
        Returns:
            np.array: Inertia matrix of the rocket
        """
        return np.diag([(1/2) * total_mass * self.r_r**2,
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2),
                                       (total_mass/12) * (self.l**2 + 3*self.r_r**2)])

    
    def I_inv(self, total_mass): 
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
