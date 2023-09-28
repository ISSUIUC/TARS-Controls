import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import environment.atmosphere as atmosphere
import dynamics.motor as motor
from motor import Motor
import dynamics.forces as forces
from rocket import Rocket

class Rocket:
    motor = None
    forces = None
    stage_config = None
    
    def __init__(self, simulation_config, atm:atmosphere.Atmosphere=None, stages:list[Rocket]=None):
        self.sim_config = simulation_config
        self.cm_rocket = simulation_config["rocket"]["cm"]
        self.cm_motor = simulation_config["motor"]["cm"]

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
        
        # Add stages to rocket via this list. Only the base rocket object should have stages, each stage should be its own rocket object with no stages
        self.stages = stages
        '''
        Start with the base rocket object, which is the first stage, then increment at each separation event.
        This is used to determine which stage the rocket is currently in. Start at -1, then increment to 0, 1, etc., so that
        they align with the indices of the stages list
        '''
        self.current_stage = -1
        '''
        This is used to store the time of the latest separation event. Start at 0, then update at each separation event so that timestamp is
        based on the latest separation event
        '''
        self.separation_timestamp = 0
        
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

    def set_motor_mass(self, timestamp) -> None:
        """Sets the mass of the motor at a given time
        
        Args:
            timestamp (float): Time in seconds
        """
        if self.current_stage == -1:
            self.motor_mass = self.motor.get_mass(timestamp)
        else:
            self.motor_mass = self.stages[self.current_stage].get_motor().get_mass(timestamp - self.separation_timestamp)
        self.rocket_total_mass = self.rocket_dry_mass + self.motor_mass
    
    def get_motor_mass(self, timestamp) -> float:
        """Returns the mass of the motor at a given time
            Only used for upper stages, don't use for base stage
            use set_motor_mass for base stage
        
        Args:
            timestamp (float): Time in seconds
        
        Returns:
            float: Mass of the motor
        """
        return self.motor.get_mass(timestamp - self.separation_timestamp)
    
    def separate_stage(self, timestamp) -> None:
        """Separates the current stage from the rocket and returns the next stage
        """
        self.current_stage += 1
        self.separation_timestamp = timestamp
    
    def get_motor(self) -> Motor:
        """Returns the motor of the rocket
        """
        return self.motor
    
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

    

    def I_inv(self, total_mass) -> np.ndarray: 
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
