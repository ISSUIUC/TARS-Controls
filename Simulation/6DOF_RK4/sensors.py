import numpy as np
import util.vectors as vct
import random
import properties as prop
import atmosphere as atm

def get_accelerometer_data(x_state):
    body_accel = vct.world_to_body(*x_state[3], x_state[2])
    kx134_rms = prop.High_G_RMS * 9.8/1000 # Convert to m/(s^2)
    accel_reading = np.array([random.gauss(body_accel[0], kx134_rms), 
                      random.gauss(body_accel[1], kx134_rms), 
                      random.gauss(body_accel[2], kx134_rms)])
    return accel_reading
    
def get_gyro_data(x_state):
    body_accel = vct.world_to_body(*x_state[3], x_state[4])
    gyro_rms = prop.Gyro_RMS * np.pi / 180000
    gyro_reading = np.array([random.gauss(body_accel[0], gyro_rms), 
                             random.gauss(body_accel[1], gyro_rms), 
                             random.gauss(body_accel[2], gyro_rms)])
    return gyro_reading

def get_barometer_data(x_state):
    a = atm.Atmosphere()
    baro_rms = prop.Barometer_RMS * 100 # Pascals
    
    pressure = a.get_pressure(x_state[0, 0])
    pressure = random.gauss(pressure, baro_rms)
    altitude = a.get_altitude(pressure)
    
    return altitude
    
def get_bno_orientation(x_state):
    true_orietnation = x_state[3]
    
    bno_reading = np.array([random.gauss(true_orietnation[0], prop.Bno_error), 
                             random.gauss(true_orietnation[1], prop.Bno_error), 
                             random.gauss(true_orietnation[2], prop.Bno_error)])
    return bno_reading
    