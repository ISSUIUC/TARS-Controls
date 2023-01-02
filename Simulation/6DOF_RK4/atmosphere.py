import numpy as np

class Atmosphere:
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
        '''
        Temperature getter function based on altitude
        
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
        
    def get_pressure(self, altitude)->float:
        '''
        Pressure getter function based on altitude
        
        Args:
            altitude (float): Altitude above sea level, in meters
            
        Returns:
            pressure (float): Pressure at altitude
        '''
        altitude_h = self.get_geometric_to_geopotential(altitude) / 1000.0
        altitude_z = altitude / 1000.0
        if (altitude_h < 11) :
            pressure = 101325.0 * pow((288.15 / (288.15 - 6.5 * altitude_h)),
                                  (34.1632 / -6.5))
        

        elif (altitude_h < 20) :
            pressure = 22632.06 * np.exp(-34.1632 * (altitude_h - 11) / 216.65)
        

        elif (altitude_h < 32) :
            pressure = 5474.889 * pow((216.65 / (216.65 + (altitude_h - 20))), 34.1632)
        

        elif (altitude_h < 47) :
            pressure = 868.0187 * pow((228.65 / (228.65 + 2.8 * (altitude_h - 32))),
                                    (34.1632 / 2.8))
        

        elif (altitude_h < 51) :
            pressure = 110.9063 * np.exp(-34.1632 * (altitude_h - 47) / 270.65)
        

        elif (altitude_h < 71) :
            pressure = 66.93887 * pow((270.65 / (270.65 - 2.8 * (altitude_h - 51))),
                                    (34.1632 / -2.8))
        

        elif (altitude_h < 84.852) :
            pressure = 3.956420 * pow((214.65 / (214.65 - 2 * (altitude_h - 71))),
                                    (34.1632 / -2))
        

        elif (altitude_h < 91) :
            pressure = np.exp(0.000000 * pow(altitude_h, 4) +
                        2.159582E-06 * pow(altitude_h, 3) +
                        -4.836957E-04 * pow(altitude_h, 2) +
                        -0.1425192 * altitude_h + 13.47530)
        

        elif (altitude_z < 100) :
            pressure = np.exp(0.000000 * pow(altitude_z, 4) +
                        3.304895E-05 * pow(altitude_z, 3) +
                        -0.009062730 * pow(altitude_z, 2) +
                        0.6516698 * altitude_z + -11.03037)
        

        elif (altitude_z < 110) :
            pressure = np.exp(0.000000 * pow(altitude_z, 4) +
                        6.693926E-05 * pow(altitude_z, 3) +
                        -0.01945388 * pow(altitude_z, 2) +
                        1.719080 * altitude_z + -47.75030)
        

        elif (altitude_z < 120) :
            pressure = np.exp(0.000000 * pow(altitude_z, 4) +
                        -6.539316E-05 * pow(altitude_z, 3) +
                        0.02485568 * pow(altitude_z, 2) +
                        -3.223620 * altitude_z + 135.9355)
        

        elif (altitude_z < 150) :
            pressure = np.exp(2.283506E-07 * pow(altitude_z, 4) +
                        -1.343221E-04 * pow(altitude_z, 3) +
                        0.02999016 * pow(altitude_z, 2) +
                        -3.055446 * altitude_z + 113.5764)
        

        elif (altitude_z < 200) :
            pressure = np.exp(1.209434E-08 * pow(altitude_z, 4) +
                        -9.692458E-06 * pow(altitude_z, 3) +
                        0.003002041 * pow(altitude_z, 2) +
                        -0.4523015 * altitude_z + 19.19151)
        

        elif (altitude_z < 300) :
            pressure = np.exp(8.113942E-10 * pow(altitude_z, 4) +
                        -9.822568E-07 * pow(altitude_z, 3) +
                        4.687616E-04 * pow(altitude_z, 2) +
                        -0.1231710 * altitude_z + 3.067409)
        

        elif (altitude_z < 500) :
            pressure = np.exp(9.814674E-11 * pow(altitude_z, 4) +
                        -1.654439E-07 * pow(altitude_z, 3) +
                        1.148115E-04 * pow(altitude_z, 2) +
                        -0.05431334 * altitude_z + -2.011365)
        

        elif (altitude_z < 750) :
            pressure = np.exp(-7.835161E-11 * pow(altitude_z, 4) +
                        1.964589E-07 * pow(altitude_z, 3) +
                        -1.657213E-04 * pow(altitude_z, 2) +
                        0.04305869 * altitude_z + -14.77132)
        

        elif (altitude_z < 1000) :
            pressure = np.exp(2.813255E-11 * pow(altitude_z, 4) +
                        -1.120689E-07 * pow(altitude_z, 3) +
                        1.695568E-04 * pow(altitude_z, 2) +
                        -0.1188941 * altitude_z + 14.56718)
        

        else :
            print("Exceeding calculatable altitude!")
            pressure = -1.0
        
        # 86k to 1000k formula not sure yet
        return pressure
    
    def get_density(self, altitude)->float:
        '''
        Returns the density at a given altitude
        
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
        

        return density
    
    def get_speed_of_sound(self, altitude)->float:
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
        '''
        Returns the wind vector at a given time
        
        Args:
            tStamp (float): Time stamp in seconds
        
        Returns:
            wind_vector (float): Wind vector at time
        '''
        dir_alpha = 0.9997
        mag_alpha = 0.99
        generated_direction_variance = np.zeros(3)
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
            
        current_wind_direction_ = current_wind_direction_ / np.linalg.norm(current_wind_direction_)
        
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