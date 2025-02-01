import os
import sys


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import estimation.ekf as kf
import pandas as pandas
import matplotlib.pyplot as plt
import properties.properties as prop
import dynamics.sensors as sensors
import environment.atmosphere as Atmosphere
import util.vectors as vct
import numpy as np


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


def readData() :
   global df_lowG_timestamp, df_highG_data, df_barometer_data


   df = pandas.read_csv("Simulation/6DOF_RK4/LookUp/ekf_tester_data.csv" ,
       usecols= ['timestamp', 'highg.ax', 'highg.ay', 'highg.az', 'barometer.altitude', 'kalman.position.px', 'orientation.yaw', 'orientation.pitch', 'orientation.roll'])
   print(df[df["orientation.roll"] != None ])
   return df.to_dict()




def implementKF(measuredDict):
   global kalman_filter
   # convert dict_values objects to lists so they're easier to work with
   barometer_data = list(measuredDict["barometer.altitude"].values())
   accel_x = list(measuredDict['highg.ax'].values())
   accel_y = list(measuredDict['highg.ay'].values())
   accel_z = list(measuredDict['highg.az'].values())


   bno_ang_pos_yaw = list(measuredDict['orientation.yaw'].values)
   bno_ang_pos_pitch = list(measuredDict['orientation.pitch'].values)
   bno_ang_pos_roll = list(measuredDict['orientation.roll'].values)
  
   accel = [accel_x, accel_y, accel_z]
   bno_ang_pos = [bno_ang_pos_yaw, bno_ang_pos_pitch, bno_ang_pos_roll]
  
   # constants for testing purposes, will be updated appropriately with lookup tables
   rho = 1.292 # air density
   thrust = 7224.49 # N
   mass = 27.216
   r = 0.0508
   Cn = 0.3390199698975432132655550164
   Ca = 0.206077604507223
   Cp = 106.21885107011164421915261998
   h = 6.68
  
   kalman_filter = kf.KalmanFilter(0.01, barometer_data[0], 0, accel_x[0], 0, 0, accel_y[0], 0, 0, accel_z[0])
  
  
   for x in range(len(barometer_data)):
       R = vct.body_to_world(bno_ang_pos[x][0], bno_ang_pos[x][1], bno_ang_pos[x][2], np.eye(3))
      
       kalman_filter.priori(R, thrust, mass, r, h, Cn, Ca, Cp, rho, bno_ang_pos[x], accel[x])
       kalman_filter.update(bno_ang_pos[x], barometer_data[x], accel_x[x], accel_y[x], accel_z[x])
       kalman_dict['x'].append(kalman_filter.get_state()[0])




def plotGraph(measuredDict):
   global df_lowG_timestamp, df_highG_data, df_barometer_data
  
   # convert the dict to a list of the values in each column so the data can be plotted
   df_lowG_timestamp = measuredDict['timestamp'].values()
   df_barometer_data = measuredDict['barometer.altitude'].values()
   df_flight_state_est = measuredDict['kalman.position.px'].values()


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





