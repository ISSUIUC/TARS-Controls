import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

class Motor():
    """Represents the motor of the rocket

    Args:
        impulse: Total impulse of the motor in Ns
        mass: Total mass of the motor in kg(?)
        delay: Time delay of the motor in seconds
        lookup_file: CSV file containing the thrust curve of the motor
    """
    def __init__(self, rocket_total_mass, cm, cm_rocket, cm_motor, rocket_dry_mass, 
                 impulse, mass, delay, lookup_file: str):
        
        lookup = os.path.join(os.path.dirname(__file__), lookup_file)
        self.thrust_data = pd.read_csv(lookup)
        # print(self.thrust_data["Thrust (N)"].dtype)
        self.total_impulse = impulse
        self.total_mass = mass
        self.current_mass = mass
        self.coast_time = delay
        self.alignment = np.array([0, 0])
        self.start_time = delay
        self.cur_line = 0
        self.rocket_total_mass = rocket_total_mass
        self.cm = cm
        self.cm_rocket = cm_rocket
        self.cm_motor = cm_motor
        self.rocket_dry_mass = rocket_dry_mass


   
    def ignite(self, start_time):
        """Ignites the motor at the given time

        Args:
            start_time: Time stamp of the launch
        """
        self.start_time = start_time


    def lerp_(self, x1, x2, y1, y2, x):
        """Helper function to linearly interpolate between two values

        Args:
            x1: Lower bound of the current interval
            x2: Upper bound of the current interval
            y1: Lower bound of the mapped interval
            y2: Upper bound of the mapped interval
            x: Value to be mapped
        """
        return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)


    def get_thrust(self, time_stamp: float) -> np.ndarray:
        """Gets the thrust of the motor at a given time

        Args:
            time: Time stamp of the current state
        """
        time_stamp = time_stamp - self.start_time
        temp = 0
        # only linearly interpolate if within the bounds of the thrust curve
        if time_stamp > self.thrust_data["Time (s)"].iloc[-1] or time_stamp < self.thrust_data["Time (s)"].iloc[0]:
            temp = 0
        # find two points on the thrust curve that bound the current time and linearly interpolate
        for i in range(len(self.thrust_data)-1):
            if float(self.thrust_data["Time (s)"].iloc[i]) == time_stamp:
                temp = float(self.thrust_data["Thrust (N)"].iloc[i])
            if float(self.thrust_data["Time (s)"].iloc[i]) < time_stamp and time_stamp < float(self.thrust_data["Time (s)"][i+1]):
                temp = self.lerp_(float(self.thrust_data["Time (s)"].iloc[i]), 
                                  float(self.thrust_data["Time (s)"].iloc[i+1]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i+1]), time_stamp) #* np.cos(self.alignment)
        self.set_alignment()
        theta = np.radians(self.alignment[0])
        phi = np.radians(self.alignment[1])
        # down is positive x, phi is 0 in ideal conditions
        vector = np.array([temp * np.cos(phi), temp * np.sin(phi)
                          * np.sin(theta), temp * np.sin(phi) * np.cos(theta)])
        
        return vector
        # return temp


    def set_alignment(self) -> None:
        """Sets the alignment of the motor

        Args:
            alignment: 3D vector representing the alignment of the motor
        """
        # read current alignment from csv file
        if self.cur_line < len(self.thrust_data):
            self.alignment = np.array([float(self.thrust_data["Theta"].iloc[self.cur_line]), float(self.thrust_data["Phi"].iloc[self.cur_line])])
            self.cur_line += 1


    def get_alignment(self) -> np.ndarray:
        """Gets the alignment of the motor
        """
        return self.alignment


    def get_mass(self, time_stamp: float) -> np.float64:
        """Gets the mass of the motor at a given time
        
        Args:
            time: Time stamp of the current state
        """ 
        # shift time_stamp to reflect the time since ignition
        shifted_time = time_stamp - self.start_time
        if shifted_time < self.thrust_data["Time (s)"].iloc[0]:
            return self.total_mass
        elif self.burnout(time_stamp):
            self.current_mass = 0

        else:
            self.current_mass = self.lerp_(float(self.thrust_data["Time (s)"].iloc[0]),
                                           float(self.thrust_data["Time (s)"].iloc[-1]),
                                           self.total_mass,
                                           0,
                                           shifted_time)
        self.rocket_total_mass = self.rocket_dry_mass + self.current_mass
        self.cm = (self.cm_rocket*self.rocket_dry_mass + self.cm_motor*self.current_mass) / (self.rocket_dry_mass + self.current_mass)
        return self.current_mass

    def burnout(self, time_stamp):
        return (time_stamp - self.start_time) >= self.thrust_data["Time (s)"].iloc[-1]

    def set_coast_time(self, coast_time: np.float64) -> None:
        """Sets the coast time of the motor

        Args:
            coast_time: Time to coast after ignition
        """
        self.coast_time = coast_time

    def get_burn_time(self) -> float:
        """Gets the burn time of the motor
        """
        return float(self.thrust_data["Time (s)"].iloc[-1])
