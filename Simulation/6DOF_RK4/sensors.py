import numpy as np
import util.vectors as vct
import random
import properties as prop

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