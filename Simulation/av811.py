import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import src.kalman_filter as kalman

flight = pd.read_csv('C:/Users/kdcch/Documents/Github/flight-data/20220623/flight_computer_trimmed.csv')

#Setting up the time vector of the flight
start_timestamp = flight.at[2,'timestamp_ms']
time = np.zeros(flight.shape[0] - 2)
for i in range(2, flight.shape[0] - 20):
    time[i - 2] = flight.at[i, 'timestamp_ms'] - start_timestamp
pos_est_data = []
vel_est_data = []
accel_est_data = []
actual_sensor_pos = []
for i in range(time.size):
    pos_est_data.append(flight.at[i, 'state_est_x'] - 1300)
    vel_est_data.append(flight.at[i, 'state_est_vx'])
    accel_est_data.append(flight.at[i, 'state_est_ax'])
    #actual_sensor_pos.append(flight.at[i, 'az'])

#Mapping a timestamp to a vector of state estimate
flight_data_dict = {
    "time":time,
    "pos":pos_est_data,
    "vel": vel_est_data,
    "accel": accel_est_data
}

flight_state_estimate = {
    "time":[0],
    "pos":[0],
    "vel": [0],
    "accel": [0]
}

pos_f = 0
pos_f_noisy = 0
vel_f = 0
accel_f = 0
# time, pos, vel, accel


kalman.initialize(flight.at[2, 'state_est_x'], flight.at[2, 'state_est_vx'], flight.at[2, 'state_est_ax'], 1)
#apriori
for i in range(2, len(time) + 2):
    predicted_state = kalman.F @ kalman.x_k
    flight_state_estimate["time"].append(i - 2)
    flight_state_estimate["pos"].append(kalman.x_k[0][0])
    flight_state_estimate["accel"].append(kalman.x_k[2][0])
    flight_state_estimate["vel"].append(kalman.x_k[1][0])
    
    kalman.P_k = kalman.F @ kalman.P_k @ kalman.F.T
    #priori
    kalman.priori(0)

    kalman.update(flight.at[i, 'altitude'] - 1300, flight.at[i, 'az'], 0, 0)


'''
for i in range(0, len(time)):
    for key in flight_data_dict:
        flight_data_dict.get(key).append(flight.at[i + 2, key])
'''
print(len(time), len(kalman.kalman_dict['alt']))
plt.plot(pos_est_data, label='Estimated Position Data')
plt.plot(kalman.kalman_dict['alt'],label='Kalman Data')
plt.ylabel('Position in Meters')
plt.xlabel('Time')
plt.show()
