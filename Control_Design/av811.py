from mmap import ACCESS_COPY
import pandas as pd
import numpy as np
#from . import src

#import matplotlib as plt

flight = pd.read_csv('C:/Users/kdcch/Documents/Github/TARS-Controls/Control_Design/flight_computer_trimmed.csv')

#Setting up the time vector of the flight
start_timestamp = flight.at[2,'timestamp']
time = np.zeros(flight.shape[0] - 2)
for i in range(2, flight.shape[0]):
    time[i - 2] = flight.at[i, 'timestamp'] - start_timestamp
print(len(time))
pos_f = 0
pos_f_noisy = 0
vel_f = 0
accel_f = 0

#src.kalman_filter.initialize(pos_f, vel_f, accel_f, time[0])

