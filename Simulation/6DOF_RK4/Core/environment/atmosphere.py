import numpy as np
import util.vectors as vct
from util.random_noise import Perlin

class Atmosphere:
    perlin = Perlin()
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
        
    def get_geometric_to_geopotential(self, altitude)->float:
        r = 6371000.0 # Radius of Earth
        return (r*altitude)/(r+altitude)
    
    def get_temperature(self, altitude):
        '''Temperature getter function based on altitude
        
        Args:
            altitude (float): Altitude above sea level, in meters
        
        Returns:
            temperature (float): Temperature at altitude
        '''
        altitude_h = self.get_geometric_to_geopotential(altitude) / 1000.0 # Geopotential altitude in km
        altitude_z = altitude / 1000.0 # Geometric altitude
        if (altitude_h < 11.0) :
            temperature = 288.15 - (6.5 * altitude_h)
        elif (altitude_h < 20.0) :
            temperature = 216.65
        

        elif (altitude_h < 32.0) :
            temperature = 196.65 + altitude_h
        

        elif (altitude_h < 47.0) :
            temperature = 139.05 + (2.8 * altitude_h)
        

        elif (altitude_h < 51.0) :
            temperature = 270.65
        

        elif (altitude_h < 71.0) :
            temperature = 413.45 - (2.8 * altitude_h)
        

        elif (altitude_h < 84.852) :
            temperature = 356.65 - (2.0 * altitude_h)
        

        elif (altitude_z < 91) :
            temperature = 186.8673
        

        elif (altitude_z < 110) :
            temperature = 263.1905 - 76.3232 * np.sqrt(1 - ((altitude_z - 91) / -19.9429)**2.0)
        

        elif (altitude_z < 120) :
            temperature = 240 + 12 * (altitude_z - 110)
        

        elif (altitude_z < 1000) :
            temperature = 1000 - 640 * np.exp(-0.01875 * ((altitude_z - 120) * (6356.766 + 120) /
                                            (6356.766 + altitude_z)))
        else :
            print("Exceeding calculatable altitude!")
            temperature = -1.0
        
        return temperature
        
    def get_pressure(self, altitude):
        """Pressure getter function based on altitude
        
        Args:
            altitude (float): Altitude above sea level, in meters
            
        Returns:
            pressure (float): Pressure at altitude
        """
        # temperature under standard condition (15 degrees C at sealevel) kelvin
        P_0 = 101325

        # pressure under standard condition in (Pa)
        T_0 = 288.16 

        # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
        b = 0.0065

        # gravitational constant 
        g = 9.81

        R = 287.05

        return P_0 * ((T_0 +(altitude)*b)/T_0)**(-g/(b*R))
    
    def get_density(self, altitude, noise=False, position=np.array([0,0,0]))->float:
        '''Returns the density at a given altitude
        
        Args:
            altitude (float): Altitude above sea level, in meters
        
        Returns:
            density (float): Density at altitude
        '''
        R = 287.053
        pressure = self.get_pressure(altitude)
        temperature = self.get_temperature(altitude)
        altitude = altitude / 1000.0

        if (altitude < 84.853) :
            density = pressure / (R * temperature)
        

        elif (altitude < 91) :
            density = np.exp(
                0.000000 * pow(altitude, 4) + -3.322622E-06 * pow(altitude, 3) +
                9.111460E-04 * pow(altitude, 2) + -0.2609971 * altitude + 5.944694)
        

        elif (altitude < 100) :
            density = np.exp(
                0.000000 * pow(altitude, 4) + 2.873405E-05 * pow(altitude, 3) +
                -0.008492037 * pow(altitude, 2) + 0.6541179 * altitude + -23.62010)
        

        elif (altitude < 110) :
            density = np.exp(
                -1.240774E-05 * pow(altitude, 4) + 0.005162063 * pow(altitude, 3) +
                -0.8048342 * pow(altitude, 2) + 55.55996 * altitude + -1443.338)
        

        elif (altitude < 120) :
            density = np.exp(
                0.00000 * pow(altitude, 4) + -8.854164E-05 * pow(altitude, 3) +
                0.03373254 * pow(altitude, 2) + -4.390837 * altitude + 176.5294)
        

        elif (altitude < 150) :
            density = np.exp(
                3.661771E-07 * pow(altitude, 4) + -2.154344E-04 * pow(altitude, 3) +
                0.04809214 * pow(altitude, 2) + -4.884744 * altitude + 172.3597)
        

        elif (altitude < 200) :
            density = np.exp(
                1.906032E-08 * pow(altitude, 4) + -1.527799E-05 * pow(altitude, 3) +
                0.004724294 * pow(altitude, 2) + -0.6992340 * altitude + 20.50921)
        

        elif (altitude < 300) :
            density = np.exp(1.199282E-09 * pow(altitude, 4) +
                        -1.451051E-06 * pow(altitude, 3) +
                        6.910474E-04 * pow(altitude, 2) + -0.1736220 * altitude +
                        -5.321644)
        

        elif (altitude < 500) :
            density = np.exp(1.140564E-10 * pow(altitude, 4) +
                        -2.130756E-07 * pow(altitude, 3) +
                        1.570762E-04 * pow(altitude, 2) + -0.07029296 * altitude +
                        -12.89844)
        

        elif (altitude < 750) :
            density = np.exp(8.105631E-12 * pow(altitude, 4) +
                        -2.358417E-09 * pow(altitude, 3) +
                        -2.635110E-06 * pow(altitude, 2) +
                        -0.01562608 * altitude + -20.02246)
        

        elif (altitude < 1000) :
            density = np.exp(-3.701195E-12 * pow(altitude, 4) +
                        -8.608611E-09 * pow(altitude, 3) +
                        5.118829E-05 * pow(altitude, 2) + -0.06600998 * altitude +
                        -6.137674)
        

        else :
            print("Exceeding calculatable altitude!")
            density = -1.0
        
        if noise:
            density *= 1 + 0.005*self.perlin.f(*(position/100))
        return density 
    
    def get_altitude(self, pressure):
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
        g = 9.81

        R = 287.05

        pressureRatio = pressure/P_0
        return -(T_0*((pressureRatio)**(b*R/(g)) - 1) * (pressureRatio)**(-b*R/(g)))/b
    
    def get_speed_of_sound(self, altitude)->float:
        '''Returns the speed of sound at a given altitude
        
        Args:
            altitude (float): Altitude above sea level, in meters
        
        Returns:
            speed of sound (float): Speed of sound at altitude
        '''
        gamma = 1.4 # Heat capacity ratio of air
        gas_constant = 287.05 # Gas constant of air
        
        return np.sqrt(gamma * gas_constant * self.get_temperature(altitude))
    
    
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
