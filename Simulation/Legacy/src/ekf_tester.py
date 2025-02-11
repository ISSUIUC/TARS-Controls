import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import estimation.ekf as kf
import pandas as pandas
import matplotlib.pyplot as plt
import properties.properties as prop

df_lowG_timestamp = 0
df_highG_data = 0
df_barometer_data = 0

# measurementDataDict = {
#     "ax": [],
#     "barometer_altitude": [],
#     "timestamp": []
# }
kalman_filter = None
kalman_dict = {"x": []}

def readData () :
    global df_lowG_timestamp, df_highG_data, df_barometer_data

    df = pandas.read_csv("../lookup/flight_computer_20221029.csv" ,
        usecols= ['timestamp_ms', 'highg_ax', 'highg_ay', 'highg_az', 'barometer_altitude', 'state_est_x'])

    return df.to_dict()


def implementKF (measuredDict):
    global kalman_filter
    # convert dict_values objects to lists so they're easier to work with
    barometer_data = list(measuredDict["barometer_altitude"].values())
    accel_x = list(measuredDict['highg_ax'].values())
    accel_y = list(measuredDict['highg_ay'].values())
    accel_z = list(measuredDict['highg_az'].values())
    kalman_filter = kf.KalmanFilter(0.01, barometer_data[0], 0, accel_x[0]*9.8, 0, 0, accel_y[0]*9.8, 0, 0, accel_z[0]*9.8)
    for x in range(len(barometer_data)):
        kalman_filter.priori()
        kalman_filter.update(barometer_data[x], accel_x[x], accel_y[x]*-9.8, accel_z[x]*-9.8)
        kalman_dict['x'].append(kalman_filter.get_state()[0])


def plotGraph (measuredDict):
    global df_lowG_timestamp, df_highG_data, df_barometer_data
    
    # convert the dict to a list of the values in each column so the data can be plotted
    df_lowG_timestamp = measuredDict['timestamp_ms'].values()
    df_barometer_data = measuredDict['barometer_altitude'].values()
    df_flight_state_est = measuredDict['state_est_x'].values()

    plt.plot(df_lowG_timestamp, df_barometer_data, label = "barometer data")
    plt.plot(df_lowG_timestamp, kalman_dict['x'], label = "state estimate")
    plt.plot(df_lowG_timestamp, df_flight_state_est, label = "flight state estimate")


    plt.title('State estimate vs Measurement')
    plt.xlabel('timestamp (ms)')
    plt.ylabel('altitude (m)')
    plt.legend()
    plt.show()



##### CALLING THE METHODS ##########
# Only call read data once to save some computational time
data = readData()
implementKF(data)
plotGraph(data)

