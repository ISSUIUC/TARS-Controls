import numpy as np
import plotly.graph_objects as go
import pandas as pd 
import statistics

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

    
def main():
    arr = pd.read_csv("Sensor Analysis\sensor_analysis.csv")
    arr_df = arr.values #changes to an array
    # arr_df = np.transpose(arr_df)
    dict_sensor = {
        "t":arr_df[:,0],
        #"alt": arr_df[:,9],
        #"highg_ax": arr_df[:,15],
        #"highg_ay": arr_df[:,16],
        #"highg_az": arr_df[:,17],
        "gx": arr_df[:,4],
        #"gy": arr_df[:,5],
        #"gz": arr_df[:,6],
        #"bno_ax": arr_df[:,20],
        #"bno_ay": arr_df[:,21],
        #"bno_az": arr_df[:,22],
        #"bno_gx": arr_df[:,23],
        #"bno_gy": arr_df[:,24],
        #"bno_gz": arr_df[:,25],
        #"bno_mx": arr_df[:,26],
        #"bno_my": arr_df[:,27],
        #"bno_mz": arr_df[:,28],
        #"magnet_mx": arr_df[:,32],
        #"magnet_my": arr_df[:,33],
        #"magnet_mz": arr_df[:,34]
        # "bno_pitch": [arr_df[:,29]],
        # "bno_roll": [arr_df[:,30]],
        # "bno_yaw": [arr_df[:,31]],
    }
    
    ## create a for loop to set t1 and t2 for whole value of times throughout launch
    ## find std dev for each t1 to t2 set
    ## add all std dev to an array
    ## plot array on graph
    # (2) Plot the data via a LRQ / Plot statistcal model

    for key in dict_sensor:
        stddev_plot(arr_df[:,0],dict_sensor[key])
        fig = go.Figure()
        ## Plotting of each dictionary entry
        fig.add_trace(go.Scatter(x=arr_df[:,0],y=dict_sensor[key],))
        
        ## Plotting Standard deviation
        fig.add_trace(go.Scatter(x=std_plot,y=time_plot,name = std_dev))
        # changing to float because dictionary values are in object datatype

        m,b = np.polyfit((dict_sensor["t"]).astype(float), (dict_sensor[key]).astype(float), 1)
        ## Line of best fit
        polyfit1 = [m * x + b for x in dict_sensor["t"]]
        fig.add_trace(go.Scatter(x=dict_sensor["t"], y=polyfit1, mode = 'lines', name = key))
        #red line in the center 
        # Finding r^2
        corr_matrix = np.corrcoef(dict_sensor[key].astype(float), polyfit1)
        corr = corr_matrix[0,1]
        R_sq = corr**2

        # def sensorlist():
            # take input based on the sensors that the user wants graphs of , and then work from there
    
        
        ### Standard Deviation Function- for a certain time range 
        

        
        print(R_sq)
        std_sample = statistics.stdev(dict_sensor[key])
        fig.add_annotation(
            text=f'<b> σ = {std_sample:.5f}</b>',
            x = 4100000,
            y = 4,
            showarrow=False,
            font=dict(size=16,color="black"),
            
        )
        fig.add_annotation(
             text=f'<b>R² = {R_sq:.5f}</b>',
             x=3900000,  # Adjust the x and y coordinates as needed
             y=4,
            showarrow=False,
            font=dict(size=16,color="black")
        )
        fig.show()
        ## use formula for sample standard deviation, create function to find standard deviation for the vars, subtract std dev factor times time from the data and see what there is
        

## assuming syntax for the line of best fit is identical to numpy

## for the line of best fit - a,b = np.polyfit(x,y,1)
## for the scatter plot ( where x and y are the arrays with data)-->
## plt.scatter(x, y) -- > adding points to plot 
## adding line of best fit plot -->
## plt.plot(x, a*x +  b, and then customizations based on how you want the line to look)

def stddev(t1, t2,time_a2,sensorkey2):
            # create the array from index t1 to t2 
            tarray = []
            keyarray = []
            # for t in range(t1,t2):
            #    tarray = tarray.append(time_a2)   
            keyarray = sensorkey2
            slicedkey = keyarray[t1:t2]
            return (statistics.stdev(slicedkey))
def stddev_plot(time_a,sensorkey):
    timearray = time_a
    std_plot = []
    time_plot = []
    for i in range(0,len(timearray)-10):
        t1 = i
        t2 = i+10 
        std_plot = std_plot.append(stddev(t1,t2,time_a,sensorkey))
        time_plot = time_plot.append((t1+t2)/2)
    
    return (std_plot,time_plot)


if (__name__ == "__main__"):
    main()