import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import environment.atmosphere as atmosphere
import dynamics.motor as motor
from dynamics.motor import Motor
import dynamics.forces as forces

class Rocket:
    motor = None
    forces = None
    stage_config = None
    
    def __init__(self, stage_config, atm:atmosphere.Atmosphere=None, stages:list=[]):
        self.stage_config = stage_config
        self.cm_rocket = stage_config["rocket_body"]["structure_cm"]
        self.cm_motor = stage_config["motor"]["cm"]
        self.cm = stage_config['rocket_body']['combined_cm']
        self.cp = stage_config['rocket_body']['combined_cp']

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
    
    def separate_stage(self, timestamp) -> bool:
        """ "Separates" the current stage from the rocket (increments the stage counter)
            and updates the separation timestamp
            
            Retuns True if there are more stages to separate, False if there are no more stages to separate
        """
        if self.current_stage == len(self.stages) - 1:
            return False
        self.current_stage += 1
        self.separation_timestamp = timestamp
        return True
            
    def get_motor(self) -> Motor:
        """Returns the motor of the rocket
        """
        return self.motor if self.current_stage == -1 else self.stages[self.current_stage].get_motor()
    
    def get_rocket_dry_mass(self) -> float:
        """Returns the dry mass of the rocket
        """
        return self.rocket_dry_mass if self.current_stage == -1 else self.stages[self.current_stage].get_rocket_dry_mass()
    
    def get_CM(self):
        return self.cm if self.current_stage == -1 else self.stages[self.current_stage].get_CM()
    
    def get_CP(self):
        return self.cp if self.current_stage == -1 else self.stages[self.current_stage].get_CP()
    
    def get_A(self):
        return self.A if self.current_stage == -1 else self.stages[self.current_stage].get_A()
    
    def get_A_s(self):
        return self.A_s if self.current_stage == -1 else self.stages[self.current_stage].get_A_s()
    
    def get_Rasaero(self):
        return self.rasaero if self.current_stage == -1 else self.stages[self.current_stage].get_Rasaero()
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
