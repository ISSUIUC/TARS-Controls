import numpy as np
import plotly.graph_objects as go
import pandas as pd 
import statistics
from plotly.subplots import make_subplots
## Notes:
# Hover over line of best fit to see what is being plotted without having to read the legend on the side
# Comment out the dictionary components you dont need if you want to see a particular graph closer up
# If you do this, you may need to resize the plots at the bottom of main()
# You may need to trim the flight data if you really want to use general r^2, line of best fit, or general std dev value


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

    
def main():
    arr = pd.read_csv("Sensor Analysis\sensor_analysis.csv")
    arr_df = arr.values #changes to an array
    dict_sensor = {
        "t":arr_df[:,0],
        "alt": arr_df[:,9],
        "highg_ax": arr_df[:,15],
        "highg_ay": arr_df[:,16],
        "highg_az": arr_df[:,17],
        "gx": arr_df[:,4],
        "gy": arr_df[:,5],
        "gz": arr_df[:,6],
        "bno_ax": arr_df[:,20],
        "bno_ay": arr_df[:,21],
        "bno_az": arr_df[:,22],
        "bno_gx": arr_df[:,23],
        "bno_gy": arr_df[:,24],
        "bno_gz": arr_df[:,25],
        "bno_mx": arr_df[:,26],
        "bno_my": arr_df[:,27],
        "bno_mz": arr_df[:,28],
        "magnet_mx": arr_df[:,32],
        "magnet_my": arr_df[:,33],
        "magnet_mz": arr_df[:,34],
        "bno_pitch": arr_df[:,29],
        "bno_roll": arr_df[:,30],
        "bno_yaw": arr_df[:,31],
    }
    t1 = -999
    t2 = -999
    while (isinstance(t1,int)==False or isinstance(t2,int)==False or t1==-999 or (arr_df[0,0] > t1) or (arr_df[-1,0] < t2)):
        try:
            t1 = int(input(f"What time do you want to start({arr_df[0,0]} to {arr_df[-1,0]})?: "))
        except Exception:
            pass
    while (isinstance(t1,int)==False or isinstance(t2,int)==False or t2== -999 or (arr_df[0,0] > t1) or (arr_df[-1,0] < t2)):
        try:
            t2 = int(input(f"What time do you want to stop({arr_df[0,0]} to {arr_df[-1,0]})?: "))
        except Exception:
            pass
    ## create a for loop to set t1 and t2 for whole value of times throughout launch
    ## find std dev for each t1 to t2 set
    ## add all std dev to an array
    ## plot array on graph
    # (2) Plot the data via a LRQ / Plot statistcal model
    fig = make_subplots(rows = len(dict_sensor), cols = 1)
    count = 0
    for key in dict_sensor:
        count = count + 1
        std_plot, time_plot, data_plot = stddev_plot(arr_df[:,0],dict_sensor[key],t1,t2)
        # fig = go.Figure()
        ## Plotting of each dictionary entry
        fig.add_trace(go.Scatter(x=time_plot,y=data_plot, name='raw data: ' + key),row=count,col=1)
        
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
        if (max(data_plot)/2 > 1):
            y_val = max(data_plot)/2
        else:
            y_val = 2
        fig.add_annotation(
            text=f'<b> σ = {std_sample:.5f}</b>',
            x = (time_plot[0]+time_plot[-1])/2,
            y = y_val,
            showarrow=False,
            font=dict(size=16,color="black"),
            row=count,
            col=1,
        )
        fig.add_annotation(
             text=f'<b>R² = {R_sq:.5f}</b>',
             x=(time_plot[0]+time_plot[-1])/2,  # Adjust the x and y coordinates as needed
             y=0,
            showarrow=False,
            font=dict(size=16,color="black"),
            row=count,
            col=1,
        )
        # fig.show()
        ## use formula for sample standard deviation, create function to find standard deviation for the vars, subtract std dev factor times time from the data and see what there is
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
    slicedkey = keyarray[t1:t2]
    return (statistics.stdev(slicedkey))
def stddev_plot(time_a,sensorkey,t_1,t_2):
    index1 = 0
    index2 = 0
    for j in range(0,len(time_a)):
        if (int(time_a[j]) - int(t_1) <= 100):
            index1 = j
        if (int(time_a[j]) - int(t_2) <= 100):
            index2 = j
    timearray = time_a[index1:index2]
    keyarray = sensorkey[index1:index2]
    # could probably be optimized with numpy arrays if runtime becomes a problem
    std_plot = np.zeros(len(timearray)-50)
    time_plot = np.zeros(len(timearray)-50)
    for i in range(0,len(timearray)-50):
        # choose t1 and t2 to move across the data set and find a std dev for each interval that will then be
        # plotted as a continuous function
        t1 = i
        t2 = i+50 
        std_plot[i]=(stddev(t1,t2,time_a,sensorkey))
        # time associated with std dev calculation is in between the two times used for the individual std dev calculation
        time_plot[i]=((timearray[t1]+timearray[t2])/2)
        
    ## returns array of standard deviation measurements and the times corresponding to them
    return (std_plot,time_plot,keyarray)


if (__name__ == "__main__"):
    main()