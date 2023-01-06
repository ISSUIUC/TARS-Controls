import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import properties as pr


"""
Represents the motor of the rocket

Args:
    impulse: Total impulse of the motor in Ns
    mass: Total mass of the motor in kg(?)
    delay: Time delay of the motor in seconds
    lookup_file: CSV file containing the thrust curve of the motor
"""
class Motor():
    def __init__(self, impulse=pr.impulse, mass=pr.motor_mass, delay=pr.delay, lookup_file: str=pr.lookup_file):
        lookup = os.path.join(os.path.dirname(__file__), lookup_file)
        self.thrust_data = pd.read_csv(lookup)
        # print(self.thrust_data["Thrust (N)"].dtype)
        self.total_impulse = impulse
        self.total_mass = mass
        self.current_mass = mass
        self.coast_time = delay
        self.alignment = np.array([0, 0])
        self.start_time = 0


    """
    Ignites the motor at the given time

    Args:
        start_time: Time stamp of the launch
    """
    def ignite(self, start_time):
        self.start_time = start_time


    """
    Helper function to linearly interpolate between two values

    Args:
        x1: Lower bound of the current interval
        x2: Upper bound of the current interval
        y1: Lower bound of the mapped interval
        y2: Upper bound of the mapped interval
        x: Value to be mapped
    """
    def lerp_(self, x1, x2, y1, y2, x):
        return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)


    """
    Gets the thrust of the motor at a given time

    Args:
        time: Time stamp of the current state
    """
    def get_thrust(self, time: float):
        # if (time < 7/.01):
        #     return np.array([2500.,0.,0.])
        # else:
        #     return np.array([0.,0.,0.])
        # shift time to reflect the time since ignition
        time = (time - self.start_time)*.01
        temp = 0
        # only linearly interpolate if within the bounds of the thrust curve
        if time > self.thrust_data["Time (s)"].iloc[-1] or time < self.thrust_data["Time (s)"].iloc[0]:
            temp = 0
        # find two points on the thrust curve that bound the current time and linearly interpolate
        for i in range(len(self.thrust_data)-1):
            if float(self.thrust_data["Time (s)"].iloc[i]) == time:
                temp = float(self.thrust_data["Thrust (N)"].iloc[i])
            if float(self.thrust_data["Time (s)"].iloc[i]) < time and time < float(self.thrust_data["Time (s)"][i+1]):
                temp = self.lerp_(float(self.thrust_data["Time (s)"].iloc[i]), 
                                  float(self.thrust_data["Time (s)"].iloc[i+1]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i+1]), time) #* np.cos(self.alignment)
        return np.array([temp,0,0])


    """
    Sets the alignment of the motor

    Args:
        alignment: 3D vector representing the alignment of the motor
    """
    def set_alignment(self, alignment: np.ndarray) -> None:
        self.alignment = alignment


    """
    Gets the alignment of the motor
    """
    def get_alignment(self) -> np.ndarray:
        return self.alignment


    """
    Gets the mass of the motor at a given time
    
    Args:
        time: Time stamp of the current state
    """ 
    def get_mass(self, time: float) -> np.float64:
        # shift time to reflect the time since ignition
        time = time - self.start_time
        if time < self.thrust_data["Time (s)"].iloc[0]:
            return self.total_mass
        elif time >= self.thrust_data["Time (s)"].iloc[-1]:
            self.current_mass = 0
        else:
            self.current_mass = self.lerp_(float(self.thrust_data["Time (s)"].iloc[0]),
                                           float(self.thrust_data["Time (s)"].iloc[-1]),
                                           self.total_mass,
                                           0,
                                           time)
        return self.current_mass


    """
    Sets the coast time of the motor

    Args:
        coast_time: Time to coast after ignition
    """
    def set_coast_time(self, coast_time: np.float64) -> None:
        self.coast_time = coast_time


if __name__ == '__main__':
    motor = Motor(69, 420, 21, '../Lookup/Cesaroni_17907N2540-P_Trimmed.csv')
    time = np.linspace(0, 10, 1000)
    # print(time)
    # print(motor.get_thrust(0))
    # thrust = np.array([motor.get_thrust(t) for t in time], dtype=np.float64)
    thrust = np.array([motor.get_mass(t) for t in time], dtype=np.float64)
    # print(thrust.shape)
    # print(time.shape)
    # print(thrust)
    # print(motor.lerp_(0.0, 1.0, 0.0, 1.0, 0.5))
    plt.plot(time, thrust)
    plt.show()