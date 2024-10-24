import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import estimation.ekf as kf
import pandas as pd
import matplotlib.pyplot as plt
import properties.properties as prop
import seaborn as sns

#plot y and z
#plot v

# measurementDataDict = {
#     "ax": [],
#     "barometer_altitude": [],
#     "timestamp": []
# }
kalman_filter = None
kalman_dict = {"x": []}

def readData () :
    global barometer_data, highG_data, lowG_data

    df = pd.read_csv("Simulation/6DOF_RK4/LookUp/data.csv" ,
        usecols= ['timestamp', 'highg.ax', 'highg.ay', 'highg.az', 'barometer.altitude', 'orientation.roll','orientation.pitch','orientation.yaw'])

    return df.to_dict()




predicted_altitude = []

def implementKF (measuredDict):
    global kalman_filter, predicted_altitude, barometer_data
    print(4)
    # convert dict_values objects to lists so they're easier to work with
    barometer_data = list(measuredDict["barometer.altitude"].values())
    accel_x = list(measuredDict['highg.ax'].values())
    accel_y = list(measuredDict['highg.ay'].values())
    accel_z = list(measuredDict['highg.az'].values())
    bno_attitude_roll = list(measuredDict['orientation.roll'].values())
    bno_attitude_pitch = list(measuredDict['orientation.pitch'].values())
    bno_attitude_yaw =list(measuredDict['orientation.yaw'].values())


    last_valid_altitude = 0
    last_valid_accel_x = 0
    last_valid_accel_y = 0
    last_valid_accel_z = 0
    last_valid_roll = 0
    last_valid_pitch = 0
    last_valid_yaw = 0  


    kalman_filter = kf.KalmanFilter(0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    for i in range(0,600000):
        if i%10000 == 0: print(i)

        if not pd.isna(barometer_data[i]):
            #print(altitude[i])
            last_valid_altitude = barometer_data[i]
        if not pd.isna(accel_x[i]):
            last_valid_accel_x = accel_x[i]
        if not pd.isna(accel_y[i]):
            last_valid_accel_y = accel_y[i]
        if not pd.isna(accel_z[i]):
            last_valid_accel_z = accel_z[i]
        if not pd.isna(bno_attitude_roll[i]):
            last_valid_roll = bno_attitude_roll[i]
        if not pd.isna(bno_attitude_pitch[i]):
            last_valid_pitch = bno_attitude_pitch[i]
        if not pd.isna(bno_attitude_yaw[i]):
            last_valid_yaw = bno_attitude_yaw[i]

        kalman_filter.priori()


        kalman_filter.update((last_valid_roll,last_valid_pitch,last_valid_yaw),last_valid_altitude, last_valid_accel_x, last_valid_accel_y*-9.8, last_valid_accel_z*-9.8)
        state = kalman_filter.get_state()
        predicted_altitude.append(state[0]) 

        


def plotGraph (measuredDict):
    global barometer_data,predicted_altitude
    fig, axes = plt.subplots(1, 2, sharex=True, figsize=(12,5))
    fig.suptitle("Dgarg")
    axes[0].set_title("Predicting altitude")
    axes[1].set_title("Barometer data")
    df_pred_altitude = pd.DataFrame(predicted_altitude, columns =['altitude'])
    df_pred_barometer_data = pd.DataFrame(barometer_data, columns =['altitude'])

    sns.lineplot(ax = axes[0], x=df_pred_altitude.index, y='altitude', data=df_pred_altitude, label='Predicted Altitude (Kalman)')
    sns.lineplot(ax = axes[1], x=df_pred_barometer_data.index, y='altitude', data=df_pred_barometer_data, label='Barometer Data (Kalman)')
    # plt.plot(df_lowG_timestamp, df_flight_state_est, label = "flight state estimate")
    print(df_pred_altitude.shape)

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

