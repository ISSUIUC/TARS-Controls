import os
import sys
import pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import estimation.ekf as ekf
import pandas as pandas
import matplotlib.pyplot as plt
import properties.properties as prop
import properties.data_loader as dataloader
import dynamics.sensors as sensors
import environment.atmosphere as atmosphere
import util.vectors as vct
import numpy as np
import dynamics.rocket as rocket_model


#### PLEASE READ BELOW COMMENTS #####

# ekf_tester_data.csv is too big for github to handle it seems (can't push), so when running this program
# you have to manually upload the file when running. Change the file name in readData() accordingly.

# original ekf_tester is Simulation/Legacy/src/ekf_tester.py


df_lowG_timestamp = 0
df_highG_data = 0
df_barometer_data = 0

kf = None
kalman_dict = {"x": []}

dir = os.path.dirname(os.path.dirname(os. path.abspath(__file__)))
csv_path = os.path.join(dir, "LookUp", "ekf_cd_test.csv")


def readData() :
   global df_lowG_timestamp, df_highG_data, df_barometer_data
   
   df = pandas.read_csv("Simulation/6DOF_RK4/LookUp/ekf_tester_data.csv" ,
       usecols= ['timestamp', 'highg.ax', 'highg.ay', 'highg.az', 'barometer.altitude', 'kalman.position.px', 'orientation.yaw', 'orientation.pitch', 'orientation.roll'])
   #print(df[df["orientation.roll"] != None ])
   
   return df


# TODO: Plotting works -> format and label properly, then verification

def implementKF(measuredDict):
   global kf
   
   # convert dict_values objects to lists so they're easier to work with
   barometer_data = list(measuredDict["barometer.altitude"].dropna())
   # print("Printing baro data: ", barometer_data[:100])
   accel_x = list(measuredDict['highg.ax'].dropna())
   accel_y = list(measuredDict['highg.ay'].dropna())
   accel_z = list(measuredDict['highg.az'].dropna())

   bno_ang_pos_yaw = list(measuredDict['orientation.yaw'].dropna())
   bno_ang_pos_pitch = list(measuredDict['orientation.pitch'].dropna())
   bno_ang_pos_roll = list(measuredDict['orientation.roll'].dropna())
  
   # unchanging constants
   rho = 1.225 # air density
   r = 0.0508 # rocket radius 
   h = 6.68 # rocket height 
   
   # constants for testing purposes, will be updated appropriately with lookup tables
   #thrust = [7224.49, 0, 0] # N #incorrect 
   # mass = 27.216 #incorrect 
   # Cn = 0.3390199698975432132655550164 #incorrect 
   # Ca = 0.206077604507223 #incorrect 
   # Cp = 106.21885107011164421915261998 #incorrect 
   
   x0 = np.zeros((6, 3))
   x0[3] = [0, 0.05, 0]
   dt = 0.01
   atm = atmosphere.Atmosphere(enable_direction_variance=True, enable_magnitude_variance=True)

   stages = []
   for stage in config['rocket']['stages'][1:]:
      stages.append(rocket_model.Rocket(dt, x0, stage, atm=atm))
   print(stages)
   rocket = rocket_model.Rocket(dt, x0, config['rocket']['stages'][0], atm=atm, stages=stages)
   
   rocket.coeffs_df = pd.read_csv(csv_path)
   rocket.coeffs_dict = {
            "CN": [0],
            "CA Power-On": [0],
            "CA Power-Off": [0],
            "CD Power-On": [0],
            "CD Power-Off": [0],
            "CL": [0],
            "CP": [0]
        }

   kf = ekf.KalmanFilter(0.01, barometer_data[0], 0, accel_x[0], 0, 0, accel_y[0], 0, 0, accel_z[0])
  
   for x in range(len(barometer_data)):
      timestamp = x * dt
      
      thrust = rocket.get_motor().get_thrust(timestamp)
      mass = rocket.get_motor_mass(timestamp)
      
      state = kf.get_state()
      
      velocity_arr = [state[1], state[4], state[7]] 
        
      if not any(component is None for component in velocity_arr):
         velocity = np.linalg.norm([state[1], state[4], state[7]])
         print(f"Velocity: {velocity}")

         try:
            rocket.update_coeffs(velocity)
         except IndexError:
            print(f"Invalid Mach number for velocity: {velocity}. Skipping coefficient update.")
            continue
      else:
         print("One or more velocity components are None. Skipping coefficient update.")
         continue
      
      Cn = rocket.get_cn()
      Ca = rocket.get_ca_on()
      Cp = rocket.get_cp()
      
      R = vct.body_to_world(bno_ang_pos_yaw[x], bno_ang_pos_pitch[x], bno_ang_pos_roll[x], np.eye(3))
      
      bno_ang_pos = (bno_ang_pos_yaw[x], bno_ang_pos_pitch[x], bno_ang_pos_roll[x])
      accel = (accel_x[x], accel_y[x], accel_z[x])
      
      kf.priori(R, thrust, mass, r, h, Cn, Ca, Cp, rho, bno_ang_pos, accel)
      kf.update(bno_ang_pos, barometer_data[x], accel_x[x], accel_y[x], accel_z[x])

      kalman_dict['x'].append(kf.get_state()[0])


def plotGraph(measuredDict):
   global df_lowG_timestamp, df_highG_data, df_barometer_data
  
   # convert the dict to a list of the values in each column so the data can be plotted
   df_lowG_timestamp = list(measuredDict['timestamp'].dropna())
   df_barometer_data = list(measuredDict['barometer.altitude'].dropna())
   df_flight_state_est = list(measuredDict['kalman.position.px'].dropna())
   df_accel_x = list(measuredDict['highg.ax'].dropna())
   df_accel_y = list(measuredDict['highg.ay'].dropna())
   df_accel_z = list(measuredDict['highg.az'].dropna())


   
   fig_kalman, (state_est, barometer, accel_x) = plt.subplots(3,1,figsize=(15,10), sharex=True)
   
   fig_kalman.suptitle('Kalman Filter Position, Barometer Altitude, and Acceleration')
   plt.xlabel("Time (s)", fontsize = 14)
      

   barometer.plot(df_lowG_timestamp[:10000], df_barometer_data[:10000], label = "Barometer Data")
   state_est.plot(df_lowG_timestamp[:6368], kalman_dict['x'][:6368], label = "State Estimate")
   accel_x.plot(df_lowG_timestamp[:1179038], df_accel_x[:1179038], label = "Acceleration (x)")
   # plt.plot(df_lowG_timestamp, df_flight_state_est, label = "Flight State Estimate")
   # plt.tight_layout()
   plt.show()
   
   # plt.title('Accel vs Time')
   # plt.xlabel('Timestamp (ms)')
   # plt.ylabel('Accel (m/s)')
   # plt.legend()
   # plt.show()


# Only call read data once to save some computational time

config = dataloader.config

if __name__ == '__main__':
  
   data = readData()
   implementKF(data)
   plotGraph(data)

