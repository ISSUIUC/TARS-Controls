import numpy as np
import util.vectors as vct
import random
import properties.properties as prop
import properties.data_loader as dataloader
import environment.atmosphere as atm

def get_accelerometer_data(x_state, sensor_config):
    """Returns the accelerometer data in the body frame of the rocket.
    
    Args:
        x_state: The state of the rocket in the world frame.
    
    Returns:
        accel_reading (np.ndarray): The accelerometer data in the body frame of the rocket. (Specific force)
    """
    gravity = [(prop.G*prop.m_e)/((prop.r_e+x_state[0,0])**2), 0, 0]
    body_accel = vct.world_to_body(*x_state[3], x_state[2] + gravity)
    kx134_rms = sensor_config["high_g"]["RMS"] * 9.81/1000 # Convert to m/(s^2)
    accel_reading = np.array([random.gauss(body_accel[0], kx134_rms), 
                      random.gauss(body_accel[1], kx134_rms), 
                      random.gauss(body_accel[2], kx134_rms)])
    return accel_reading
    
def get_gyro_data(x_state, sensor_config):
    """Returns the gyroscope data in the body frame of the rocket.
    
    Args:
        x_state: The state of the rocket in the world frame.
        
    Returns:
        gyro_reading (np.ndarray): The gyroscope data in the body frame of the rocket. (Angular velocity)
    """
    body_accel = vct.world_to_body(*x_state[3], x_state[4])
    gyro_rms = sensor_config["gyro"]["RMS"] * np.pi / 180000
    gyro_reading = np.array([random.gauss(body_accel[0], gyro_rms), 
                             random.gauss(body_accel[1], gyro_rms), 
                             random.gauss(body_accel[2], gyro_rms)])
    return gyro_reading

def get_barometer_data(x_state, sensor_config):
    """Returns the barometer data in the world frame of the rocket.
    
    Args:
        x_state: The state of the rocket in the world frame.
    
    Returns:
        altitude (float): The altitude of the rocket in meters. (Barometric)
    """
    a = atm.Atmosphere()
    baro_rms = sensor_config["barometer"]["RMS"] * 100 # Pascals
    
    pressure = a.get_pressure(x_state[0, 0])
    pressure = random.gauss(pressure, baro_rms)
    altitude = a.get_altitude(pressure)
    
    return altitude
    
def get_bno_orientation(x_state, sensor_config):
    """Returns the sensor emulated orientation data in the world frame of the rocket.
    
    Args:
        x_state: The state of the rocket in the world frame.
        
    Returns:
        bno_reading (np.ndarray): The sensor emulated orientation data in the world frame of the rocket. (Euler angles)
    """
    true_orientation = x_state[3]
    bno_error = sensor_config["bno"]["error"] * np.pi / 180
    bno_reading = np.array([random.gauss(true_orientation[0], bno_error), 
                             random.gauss(true_orientation[1], bno_error), 
                             random.gauss(true_orientation[2], bno_error)])
    return bno_reading
    