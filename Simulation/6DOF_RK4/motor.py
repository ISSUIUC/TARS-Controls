import numpy as np
import pandas as pd

class Motor():
    def __init__(self, impulse, mass, delay, lookup_file: str):
        lookup = os.path.join(os.path.dirname(__file__), lookup_file)
        self.thrust_data = pd.read_csv(lookup)
        print(self.thrust_data["Thrust (N)"].dtype)
        self.total_impulse = impulse
        self.current_mass = mass
        self.coast_time = delay
        self.alignment = np.array([0, 0])
  
    def lerp_(self, x1, x2, y1, y2, x):
        return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)

    def get_thrust(self, time: float):
        # given a CSV, find which two timesteps time is between and interpolate between the two timestamps to find the thrust
        if time > self.thrust_data["Time (s)"].iloc[-1] or time < self.thrust_data["Time (s)"].iloc[0]:
            return 0
        for i in range(len(self.thrust_data)-1):
            if float(self.thrust_data["Time (s)"].iloc[i]) == time:
                return float(self.thrust_data["Thrust (N)"].iloc[i])
            if float(self.thrust_data["Time (s)"].iloc[i]) < time and time < float(self.thrust_data["Time (s)"][i+1]):
                return self.lerp_(float(self.thrust_data["Time (s)"].iloc[i]), 
                                  float(self.thrust_data["Time (s)"].iloc[i+1]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i]), 
                                  float(self.thrust_data["Thrust (N)"].iloc[i+1]), time) #* np.cos(self.alignment)
  
    def set_alignment(self, alignment: np.ndarray) -> None:
        self.alignment = alignment
        
    def get_alignment(self) -> np.ndarray:
        return self.alignment
        
    # def get_mass() -> np.float64:
        
    # def set_coast_time(coast_time: np.float64) -> None:
    #     coast_time = 

import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    motor = Motor(69, 420, 21, '..\Lookup\Cesaroni_17907N2540-P_Trimmed.csv')
    time = np.linspace(0, 10, 1000)
    # print(time)
    # print(motor.get_thrust(0))
    thrust = np.array([motor.get_thrust(t) for t in time], dtype=np.float64)
    # print(thrust.shape)
    # print(time.shape)
    # print(thrust)
    # print(motor.lerp_(0.0, 1.0, 0.0, 1.0, 0.5))
    plt.plot(time, thrust)
    plt.show()