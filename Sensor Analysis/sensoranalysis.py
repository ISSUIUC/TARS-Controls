import numpy as np
import plotly.graph_objects as go
# import analysis.csv as csv
# --> steps to perform for the analysis
## --> code corresponding to steps
# Altitude (potentially), highg_ax, highg_ay, highg_az, gx, gy, gz, bno_ax, bno_ay, bno_az, bno_gx, bno_gy, bno_gz, 
# bno_mx, bno_my, bno_mz, magnet_mx, magnet_my, magnet_mz, bno_pitch (potentially), bno_roll (potentially), 
# bno_yaw (potentially)

# (1) Import data and create a dictionary for the flight data 
# dict = {
#     "t":[];
#     "alt": [];
#     "highg_ax": [];
#     "highg_ay": [];
#     "highg_az": [];
#     "gx": [];
#     "gy": [];
#     "gz": [];
#     "bno_ax": [];
#     "bno_ay": [];
#     "bno_az": [];
#     "bno_gx": [];
#     "bno_gy": [];
#     "bno_gz": [];
#     "bno_mx": [];
#     "bno_my": [];
#     "bno_mz": [];
#     "magnet_mx": [];
#     "magnet_my": [];
#     "magnet_mz": [];
#     # "bno_pitch": [];
#     # "bno_roll": [];
#     # "bno_yaw": [];
# }

# (2) Plot the data via a LRQ / Plot statistcal model?

## assuming syntax for the line of best fit is identical to numpy

## for the line of best fit - a,b = np.polyfit(x,y,1)
## for the scatter plot ( where x and y are the arrays with data)-->
## plt.scatter(x, y) -- > adding points to plot 
## adding line of best fit plot -->
## plt.plot(x, a*x +  b, and then customizations based on how you want the line to look)

