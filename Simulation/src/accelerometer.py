
import numpy.random as random
import src.constants as constants

#  accel is in m/s^2
def accelerometer_noise(accel):
    accel_g = accel / constants.g

    # LSM9DS1 IMU has an accuracy of +- 90 mG
    # 90% of test data falls within this range - calculate standard deviation associated with 90% (z = 1.645)
    # limit at ±16g, weighted switch over to high G IMU starting at 12G
    low_g_range = 90/1000 #G
    low_g_max = 16
    z = 1.645
    low_g_noise = (accel_g + random.normal(scale = low_g_range/z))* constants.g

    # KX134 IMU has an accuracy of +- (75 mG + 5mg/ºC above 25ºC), assume nosecone internal temperature is 70ºC
    # 90% of test data falls within this range - calculate standard deviation associated with 90% (z = 1.645)
    internal_temp = 70  # ºC
    room_temp = 25  # ºC
    high_g_room_temp_acc = 75 # mG
    # mg /ºC from room temperature
    conversionRate = .5

    high_g_range = (high_g_room_temp_acc + conversionRate * (internal_temp - room_temp))/1000
    high_g_noise = (accel_g + random.normal(scale = high_g_range/z))* constants.g

    crossover_range_size  = 4
    if(high_g_noise / constants.g <= low_g_max - crossover_range_size):
        return low_g_noise
    elif(high_g_noise / constants.g >= low_g_max):
        return high_g_noise
    else:
        gamma = (high_g_noise / constants.g - (low_g_max - crossover_range_size)) / crossover_range_size
        return low_g_noise * (1 - gamma) + gamma * high_g_noise