import numpy as np
import plotly.graph_objects as go
import pandas as pd 

data = pd.read_csv('/Sensor Analysis/sensor_analysis.csv') ## --> analysis csv needs to be updated 

# --> steps to perform for the analysis
## --> code corresponding to steps
# Altitude (potentially), highg_ax, highg_ay, highg_az, gx, gy, gz, bno_ax, bno_ay, bno_az, bno_gx, bno_gy, bno_gz, 
# bno_mx, bno_my, bno_mz, magnet_mx, magnet_my, magnet_mz, bno_pitch (potentially), bno_roll (potentially), 
# bno_yaw (potentially)



# (1) Import data and create a dictionary for the flight data 
#

def main():
    arr = read_csv('sensor_analysis.csv')
    arr_df = arr.values
    arr_df = np.transpose(arr_df)
    dict = {
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
        "magnet_mz": arr_df[:,34]
        # "bno_pitch": [arr_df[:,29]],
        # "bno_roll": [arr_df[:,30]],
        # "bno_yaw": [arr_df[:,31]],
    }

    # (2) Plot the data via a LRQ / Plot statistcal model

    for key in dict:
        fig = go.Figure()
        ## Plotting of each dictionary entry
        fig.add_trace(data=[go.Scatter(x=dict[key],y=arr_df[:,0],)])
        m, b = np.polyfit(dict["t"], dict[key], 1)
        ## Line of best fit
        polyfit = [m * x + b for x in dict["t"]]
        fig.add_trace(go.Scatter(dict["t"], polyfit, mode = 'lines', name = 'polyfit'))
        ## Finding r^2
        corr_matrix = np.corrcoef(dict[key], polyfit)
        corr = corr_matrix[0,1]
        R_sq = corr**2
    
        print(R_sq)
        
        fig.add_annotation(
            text=f'R² = {R_sq:.2f}',
            x=1,  # Adjust the x and y coordinates as needed
            y=4,
            showarrow=False,
            font=dict(size=12),
        )
        fig.show()
        

## assuming syntax for the line of best fit is identical to numpy

## for the line of best fit - a,b = np.polyfit(x,y,1)
## for the scatter plot ( where x and y are the arrays with data)-->
## plt.scatter(x, y) -- > adding points to plot 
## adding line of best fit plot -->
## plt.plot(x, a*x +  b, and then customizations based on how you want the line to look)



if (__name__ == "__main__"):
    main()