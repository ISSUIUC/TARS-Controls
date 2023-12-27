import numpy as np
import scipy.constants
import util.vectors as vct
from ambiance import Atmosphere as ICAOmodel
from util.random_noise import Perlin

class WindModel:
    def __init__(self, wind_direction_variance_mean = 0.0, 
                 wind_direction_variance_stddev = 0.01,
                 wind_magnitude_variance_mean = 0.0, 
                 wind_magnitude_variance_stddev = 0.5,
                 enable_direction_variance = False, 
                 enable_magnitude_variance = False,
                 nominal_wind_direction = np.array([-1.0, 0.0, 0.0]),
                 nominal_wind_magnitude = 0.0):
        
        # Wind variance generation
        self.wind_direction_variance_mean_ = wind_direction_variance_mean
        self.wind_direction_variance_stddev_ = wind_direction_variance_stddev
        self.wind_magnitude_variance_mean_ = wind_magnitude_variance_mean
        self.wind_magnitude_variance_stddev_ = wind_magnitude_variance_stddev
        
        self.enable_direction_variance_ = enable_direction_variance
        self.enable_magnitude_variance_ = enable_magnitude_variance
        
        # Nominal wind without variance
        self.nominal_wind_direction_ = nominal_wind_direction
        self.nominal_wind_magnitude_ = nominal_wind_magnitude
        
        # Variance update rate
        self.last_direction_variance_update_ = -100.0
        self.last_magnitude_variance_update_ = -100.0
        self.direction_variance_update_rate_ = 1.0
        self.magnitude_variance_update_rate_ = 1.0
        
        # Smoothed wind variance applid onto nominal wind
        self.direction_variance_vect_ = np.array([0.0, 0.0, 0.0])
        self.magnitude_variance_val_ = 0.0

    ########## Get Parameters ############
    def get_nominal_wind_direction(self) -> float:
        return self.nominal_wind_direction_
    
    def get_nominal_wind_magnitude(self) -> float:
        return self.nominal_wind_magnitude_
    
    
    ########## Set parameters ###########
    def set_nominal_wind_direction(self, direction) -> None:
        self.nominal_wind_direction_ = direction / np.linalg.norm(direction)
    
    def set_nominal_wind_magnitude(self, magnitude) -> None:
        self.nominal_wind_magnitude_ = magnitude
    
    ########## Wind modeling #################
    def toggle_wind_direction_variance(self, toggle) -> None:
        self.enable_direction_variance_ = toggle
    
    def toggle_wind_magnitude_variance(self, toggle) -> None:
        self.enable_magnitude_variance_ = toggle
        
    def get_wind_vector(self, tStamp)->np.ndarray:
        '''Returns the wind vector at a given time
        
        Args:
            tStamp (float): Time stamp in seconds
        
        Returns:
            wind_vector (float): Wind vector at time
        '''
        dir_alpha = 0.9997
        mag_alpha = 0.99
        generated_direction_variance = np.zeros(3)
        generated_magnitude_variance = 0
        current_wind_direction_ = np.zeros(3)
        if self.enable_direction_variance_:
            if((tStamp - self.last_direction_variance_update_) >= self.direction_variance_update_rate_):
                generated_direction_variance = np.random.normal(self.nominal_wind_direction_, 3)
                generated_direction_variance = generated_direction_variance / np.linalg.norm(generated_direction_variance)
                self.last_direction_variance_update_ = tStamp
            
             # Apply low-pass filter to smooth wind direction change
            self.direction_variance_vect_ = dir_alpha * self.direction_variance_vect_ + (1.0 - dir_alpha) * generated_direction_variance

            current_wind_direction = self.nominal_wind_direction_ + self.direction_variance_vect_
        else:
            current_wind_direction = self.nominal_wind_direction_
        
        current_wind_direction_ = vct.norm(current_wind_direction_)
        
        if self.enable_magnitude_variance_:
            if (tStamp - self.last_magnitude_variance_update_) >= self.magnitude_variance_update_rate_ :
                generated_magnitude_variance = np.random.normal(self.nominal_wind_magnitude_, 1)
                self.last_magnitude_variance_update_ = tStamp
            # Apply low-pass filter to smooth wind magnitude change
            self.magnitude_variance_val_ = mag_alpha * self.magnitude_variance_val_ + (1.0 - mag_alpha) * generated_magnitude_variance
            
            current_wind_magnitude = self.nominal_wind_magnitude_ + self.magnitude_variance_val_
        else:
            current_wind_magnitude = self.nominal_wind_magnitude_
            
        return current_wind_direction * current_wind_magnitude

class AtmosphereModel:
    perlin = Perlin()
    def get_temperature(altitude: float) -> float:
        """Returns the temperature at a given altitude in meters"""
        return ICAOmodel(altitude).temperature[0]
    
    def get_pressure(altitude: float) -> float:
        """Returns the pressure at a given altitude in meters"""
        return ICAOmodel(altitude).pressure[0]
    
    def get_density(altitude: float, noise=False, position=np.array([0,0,0])) -> float:
        """Returns the pressure at a given altitude in meters"""
        density = ICAOmodel(altitude).density[0]
        if noise:
            density *= 1 + 0.005*self.perlin.f(*(position/100))
        return density
    
    def get_speed_of_sound(altitude: float) -> float:
        """Returns the speed of sound at a given altitude in meters"""
        return ICAOmodel(altitude).speed_of_sound[0]
    
    def pressure_to_altitude(pressure: float) -> float:
        '''Returns altitude at a given pressure using the international barometric formula
        
        Args:
            pressure: Atmospheric pressure (pascals)
        
        Returns:
            Altitude from sea level in meters
        '''
        # temperature under standard condition (15 degrees C at sealevel) kelvin
        P_0 = 101325
        # pressure under standard condition in (Pa)
        T_0 = 288.16 
        # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
        b = 0.0065

        # gravitational constant 
        g = scipy.constants.g
        R = 287.05

        pressureRatio = pressure/P_0
        return -(T_0*((pressureRatio)**(b*R/(g)) - 1) * (pressureRatio)**(-b*R/(g)))/b

class Atmosphere:
    def __init__(self, wind_direction_variance_mean = 0.0, 
                 wind_direction_variance_stddev = 0.01,
                 wind_magnitude_variance_mean = 0.0, 
                 wind_magnitude_variance_stddev = 0.5,
                 enable_direction_variance = False, 
                 enable_magnitude_variance = False,
                 nominal_wind_direction = np.array([-1.0, 0.0, 0.0]),
                 nominal_wind_magnitude = 0.0):
        
        self.wind = WindModel(wind_direction_variance_mean, wind_direction_variance_stddev,
                                     wind_magnitude_variance_mean, wind_magnitude_variance_stddev,
                                     enable_direction_variance, enable_magnitude_variance,
                                     nominal_wind_direction, nominal_wind_magnitude)
        
