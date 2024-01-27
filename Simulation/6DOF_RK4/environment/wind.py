import numpy as np
from numpy import ndarray
import math

import atmosphere
class WindModel(atmosphere.AbstractWindModel):
    def __init__(self):
        self.altitude = []
        self.wind_speed = []
        self.wind_angle = []
        # Update the wind velocity
        # The tuple is made of 2 values, a wind velocity, and an angle in degrees for the theta
    def get_wind_vector(self, altitude, time_stamp) -> ndarray:
        # Get from the map
        velocity = np.interp(altitude, self.altitude, self.wind_speed)
        theta = np.interp(altitude, self.altitude, self.wind_angle)
        return np.array([0, math.cos(theta) * velocity, math.sin(theta) * velocity])
