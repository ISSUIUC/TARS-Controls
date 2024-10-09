import numpy as np
import plotly.graph_objects as go
import pandas as pd 
import statistics
from plotly.subplots import make_subplots
## Notes:
# Hover over line of best fit to see what is being plotted without having to read the legend on the side
# Click in the legend on the side to find a specific plot
# Comment out the dictionary components you dont need if you want to see a particular graph closer up
# If you do this, you may need to resize the plots at the bottom of main()

# stddev(t1, t2,time_a2,sensorkey2)
# program logic - 

## data = pd.read_csv('/Sensor Analysis/sensor_analysis.csv') ## --> analysis csv needs to be updated 

# --> steps to perform for the analysis
## --> code corresponding to steps
# Altitude (potentially), highg_ax, highg_ay, highg_az, gx, gy, gz, bno_ax, bno_ay, bno_az, bno_gx, bno_gy, bno_gz, 
# bno_mx, bno_my, bno_mz, magnet_mx, magnet_my, magnet_mz, bno_pitch (potentially), bno_roll (potentially), 
# bno_yaw (potentially)




# (1) Import data and create a dictionary for the flight data 

# As an amendment to this to avoid multiple tabs opening at once, create a list
## this list is defined as the sensors that you do want to use 

#also trim the flight data

""" 
    Aditya Srikanth 9-28-2024
    - Removed commented code that wasn't used with sensors or unnecessarily caused errors
    - Changed certain arr.values used in functions to dict_sensor value
        +++ arr_df[:,7] to dict_sensor["lowG_data.timeStamp_lowG"]
    - When figuring graph titles, isntead of keeping a long list of != statements, i made a list of excluded tags and checked to make sure element was not in the list
"""
    
def main():
    sensorData = pd.read_csv("Sensor Analysis/csvfinal.csv")
    arr = pd.read_csv("Sensor Analysis/csvfinal.csv")
    arr_df = np.asarray(arr.values) #changes to an array
    dict_sensor = {
        "has_lowG_data": arr_df[:,0],
        "lowG_data.ax":arr_df[:,1],
        "lowG_data.ay": arr_df[:,2],
        "lowG_data.az": arr_df[:,3],
        "lowG_data.gx": arr_df[:,4],
        "lowG_data.gy": arr_df[:,5],
        "lowG_data.gz": arr_df[:,6],
        "lowG_data.timeStamp_lowG": arr_df[:,7],
        "has_highG_data": arr_df[:,8],
        "highG_data.hg_ax": arr_df[:,9],
        "highG_data.hg_ay": arr_df[:,10],
        "highG_data.hg_az": arr_df[:,11],
        "highG_data.timeStamp_highG": arr_df[:,12],

        "has_barometer_data": arr_df[:,21],
        "barometer_data.temperature": arr_df[:,22],
        "barometer_data.pressure": arr_df[:,23],
        "barometer_data.altitude": arr_df[:,24],
        "barometer_data.timeStamp_barometer": arr_df[:,25],

        "has_voltage_data": arr_df[:,44],
        "voltage_data.v_battery": arr_df[:,45],
        "voltage_data.timestamp": arr_df[:,46],
        "has_orientation_data": arr_df[:,47],
        "orientation_data.accel.ax": arr_df[:,48],
        "orientation_data.accel.ay": arr_df[:,49],
        "orientation_data.accel.az": arr_df[:,50],
        "orientation_data.gyro.gx": arr_df[:,51],
        "orientation_data.gyro.gy": arr_df[:,52],
        "orientation_data.gyro.gz": arr_df[:,53],
        "orientation_data.magnet.mx": arr_df[:,54],
        "orientation_data.magnet.my": arr_df[:,55],
        "orientation_data.magnet.mz": arr_df[:,56],
        "orientation_data.angle.yaw": arr_df[:,57],
        "orientation_data.angle.pitch": arr_df[:,58],
        "orientation_data.angle.roll": arr_df[:,59],
        "orientation_data.timeStamp_orientation": arr_df[:,60],
        "has_magnetometer_data": arr_df[:,61],
        "magnetometer_data.magnetometer.mx": arr_df[:,62],
        "magnetometer_data.magnetometer.my": arr_df[:,63],
        "magnetometer_data.magnetometer.mz": arr_df[:,64],
        "magnetometer_data.timestamp": arr_df[:,65],
        "has_gas_data": arr_df[:,65],
        "gas_data.temp": arr_df[:,67],
        "gas_data.humidity": arr_df[:,68],
        "gas_data.pressure": arr_df[:,69],
        "gas_data.resistance": arr_df[:,70],
        "gas_data.timestamp": arr_df[:,71],

    }
    # stuff at the bottom isnt't used anywhere in this file. dont want to delete until i double check. update: made a separate commit taking this out
    
    counter = 0
    for key2 in dict_sensor.keys():
        new_array = dict_sensor[key2][dict_sensor[key2]!=0.0]
        dict_sensor[key2] = new_array
        counter += 1
    

    # preliminary time values to start the while loops
    t1 = -999
    t2 = -999
    # hello
    print("*************************************************")
    print("Please use time intervals of over 1000000        ")
    print("Data can be cropped inside if needed             ")
    print("For optimal performance, keep time intervals     ")
    print("Under the value of 10000000 (1/25th of data)     ")
    print("Average run-time (3000000 interval)  = 15 seconds")
    print("Average run-time (10000000 interval) = 75 seconds")
    print("*************************************************")
    while (isinstance(t1,int)==False or isinstance(t2,int)==False or t1==-999 or (min(sensorData["lowG_data.timeStamp_lowG"]) > t1)): #checking if time values are out of bounds or not of int type
        try:
            t1 = int(input(f'What time do you want to start({min(arr["lowG_data.timeStamp_lowG"])} to {max(sensorData["lowG_data.timeStamp_lowG"])})?: '))
        except Exception: # having the while loop run again if an integer is not entered
            pass
    while (isinstance(t1,int)==False or isinstance(t2,int)==False or t2== -999  or (max(sensorData["lowG_data.timeStamp_lowG"]) < t2)): #checking if time values are out of bounds or not of int type
        try:
            t2 = int(input(f'What time do you want to stop({min(sensorData["lowG_data.timeStamp_lowG"])} to {max(sensorData["lowG_data.timeStamp_lowG"])})?: '))
        except Exception: # having the while loop run again if an integer is not entered
            pass

    ## create a for loop to set t1 and t2 for whole value of times throughout launch
    ## find std dev for each t1 to t2 set
    ## add all std dev to an array
    ## plot array on graph
    # (2) Plot the data via a LRQ / Plot statistcal model
    titles = []
    count1 = 0
    notGraphedData = [
        "has_lowG_data", "lowG_data.timeStamp_lowG", "has_highG_data", "has_gps_data", 
        "has_barometer_data", "has_kalman_data", "has_rocketState_data", "has_flap_data", 
        "has_voltage_data", "has_orientation_data", "has_magnetometer_data", 
        "highG_data.timeStamp_highG", "gps_data.timeStamp_GPS", "barometer_data.timeStamp_barometer", 
        "kalman_data.timeStamp_state", "rocketState_data.timestamp", "flap_data.timeStamp_flaps", 
        "voltage_data.timestamp", "has_gas_data", "orientation_data.timeStamp_orientation", 
        "magnetometer_data.timestamp", "gas_data.timestamp"
    ]
    for header in dict_sensor:
        # appending only titles of things to be graphed (no time or has data columns)
        if header not in notGraphedData:
            titles.append(header)
            count1 += count1
    fig = make_subplots(rows = len(dict_sensor), cols = 1, subplot_titles=titles)
    count = 0
    count2 = 0
    for key in dict_sensor:
        count2 += 1
        # plotting only columns that contain important data
        if key not in notGraphedData:
            count += 1
            if (count2 <= 7):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,7],dict_sensor[key],t1,t2)
            if (8 <= count2 <= 12):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,12],dict_sensor[key],t1,t2)
           
            if (13 <= count2 <= 17):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,25],dict_sensor[key],t1,t2)
           
            if (18 <= count2 <= 20):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,46],dict_sensor[key],t1,t2)
            if (21 <= count2 <= 34):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,60],dict_sensor[key],t1,t2)
            if (35 <= count2 <= 39):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,65],dict_sensor[key],t1,t2)
            if (40 <= count2 <= 45):
                std_plot, time_plot, data_plot = stddev_plot(arr_df[:,71],dict_sensor[key],t1,t2)
            # fig = go.Figure()
            ## Plotting of each dictionary entry
            # Printing key to see how close to done the graphing operation is 
            print(key)
            fig.add_trace(go.Scatter(x=time_plot.astype(np.float64),y=data_plot.astype(np.float64), name='raw data: ' + key),row=count,col=1)

            ## Plotting Standard deviation
            fig.add_trace(go.Scatter(x=time_plot,y=std_plot,name = 'std_dev: ' + key),row=count,col=1)
            
            # changing to float because dictionary values are in object datatype
            m,b = np.polyfit((time_plot).astype(float), (data_plot[0:len(time_plot)]).astype(float), 1)
            
            ## Line of best fit
            polyfit1 = [m * x + b for x in time_plot]
            fig.add_trace(go.Scatter(x=time_plot, y=polyfit1, mode = 'lines', name = 'best fit: ' + key),row=count,col=1)
            
            #red line in the center 
            # Finding r^2
            corr_matrix = np.corrcoef(data_plot[0:len(time_plot)].astype(float), polyfit1)
            corr = corr_matrix[0,1]
            R_sq = corr**2

            # def sensorlist():
                # take input based on the sensors that the user wants graphs of , and then work from there
        
            
            ### Standard Deviation Function- for a certain time range 
            
            
            std_sample = statistics.stdev(data_plot)
            
            ## adding r^2 and std deviation onto the graph
            if (max(data_plot)/2 > 1):
                y_val = max(data_plot)/2
            elif (min(data_plot)/2 < 1):
                y_val = min(data_plot)/2
            else:
                y_val = 0.2
            # total standard dev over sample
            fig.add_annotation(
                text=f'<b> σ = {std_sample:.5f}</b>',
                x = (time_plot[0]+time_plot[-1])/2,
                y = y_val,
                showarrow=False,
                font=dict(size=16,color="black"),
                row=count,
                col=1,
            )
            # r^2 of best fit line over whole sample
            fig.add_annotation(
                text=f'<b>R² = {R_sq:.5f}</b>',
                x=(time_plot[0]+time_plot[-1])/2+1000000,  # Adjust the x and y coordinates as needed
                y=0,
                showarrow=False,
                font=dict(size=16,color="black"),
                row=count,
                col=1,
            )
            # fig.show()
            # use formula for sample standard deviation, create function to find standard deviation for the vars, subtract std dev factor times time from the data and see what there is
    ## change size here if you want
    fig['layout'].update(height= 5000, width=1500,
                     title='Sensor Data')
    fig.show()   

## assuming syntax for the line of best fit is identical to numpy

## for the line of best fit - a,b = np.polyfit(x,y,1)
## for the scatter plot ( where x and y are the arrays with data)-->
## plt.scatter(x, y) -- > adding points to plot 
## adding line of best fit plot -->
## plt.plot(x, a*x +  b, and then customizations based on how you want the line to look)

def stddev(t1, t2,time_a2,sensorkey2):
    # for t in range(t1,t2):
    #    tarray = tarray.append(time_a2)   
    keyarray = sensorkey2
    slicedkey = keyarray[0:t2]
    return (statistics.stdev(slicedkey))

def stddev_plot(time_a,sensorkey,t_1,t_2):
    index1 = 0
    index2 = 0
    # choosing indexes
    for j in range(0,len(time_a)):
        if (int(time_a[j]) - int(t_1) <= 20000 and time_a[j]!=0.0):
            index1 = j
        if (abs(int(time_a[j]) - int(t_2)) <= 50000 and int(time_a[j]) - int(t_2)>0):
            index2 = j
    # ensuring a final index is chosen no matter what
    if index2 == 0:
        index2 = max(time_a)
    timearray = time_a[index1:index2]
    timearray = timearray[timearray != 0.0]
    # print(len(timearray))
    keyarray = sensorkey[index1:index2]
    std_plot = np.zeros(len(timearray)-50)
    time_plot = np.zeros(len(timearray)-50)
    # print(timearray)
    for i in range(0,len(timearray)-50):
        # choose t1 and t2 to move across the data set and find a std dev for each interval that will then be
        # plotted as a continuous function
        t1 = i
        t2 = i+50 
        # std dev array of values
        std_plot[i]=(stddev(t1,t2,time_a,sensorkey))
        # time associated with std dev calculation is in between the two times used for the individual std dev calculation
        time_plot[i]=((timearray[t1]+timearray[t2])/2)
        
    ## returns array of standard deviation measurements and the times corresponding to them
    keyarray = np.array(keyarray)
    return (std_plot,time_plot,keyarray)

if (__name__ == "__main__"):
    main()
