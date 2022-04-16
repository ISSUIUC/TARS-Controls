import numpy as np
import pandas as pd
import csv
import src.atmosphere as atmosphere
import sys

# for inputs, extension is a value in inches, between 0 and 0.5. rasaero is the csv file pulled from Lookup

def cd_interpolation(altitude, velocity, ang_of_att, max_extension, extension, rasaero, before_launch, before_burnout):
    
    if (before_launch):
        return 0
    
    # Use atmosphere in src to find speed of sound given altitude
    a = atmosphere.speed_sound(altitude)

    # Define mach number for csv lookup, rounded to hundreds place
    mach = round((velocity/a), 2)
    
    # round AoA to closest integer
    ang_of_att = round(ang_of_att)

    #Define blank upper and lower Cds
    Cd_low = 0
    Cd_up = 0

    # define csv file to search through
    #csv_file = csv.reader(open('RASAero.csv', 'r'))
    csv_file = rasaero

    # Define protuberance percentage of full extension given current extension
    protub_perc = extension/max_extension

    # Define starting Cd
    Cd = 0
    
    cd_vals = csv_file["CD_Power_Off"]
    if (before_burnout):
        cd_vals = csv_file["CD_Power_On"]
    
    # Find indices where the mach values match up
    mach_indices = np.where(csv_file['mach'] == mach)[0]    
    if len(mach_indices) == 0:
        return Cd
        
    # Interpolate to find Cd value
    for idx in range(mach_indices[0], mach_indices[-1] + 1):
        if (ang_of_att == csv_file['alpha_deg'][idx] and csv_file['protuberance'][idx] <= protub_perc <= csv_file['protuberance'][idx+1]):
            
            Cd_low = cd_vals[idx]
            Cd_up = cd_vals[idx+1]

            Cd = np.interp(protub_perc, [csv_file['protuberance'][idx], csv_file['protuberance'][idx+1]], [Cd_low, Cd_up])        
    # print("Protuberance: ", protub_perc)
    # print("MACH: ", mach)
    # print("CD: ", Cd)
    # print("In Thrust: " , before_burnout)

    # print(mach_indices)

    # for loop to find upper Cd value of interpolation
    # for row in range(len(csv_file['mach'])-1):
    #     # see if row meets values for the lower Cd

    #     if (mach == csv_file['mach'][row]) & (ang_of_att == csv_file['alpha_deg'][row]) & (csv_file['protuberance'][row] <= protub_perc <= csv_file['protuberance'][row+1]):
    #         # Cd_low = csv_file['cd_power_off'][row]
    #         # Cd_up = csv_file['cd_power_off'][row+1]
            
    #         Cd_low = cd_vals[row]
    #         Cd_up = cd_vals[row+1]

    #         Cd = np.interp(protub_perc, [csv_file['protuberance'][row], csv_file['protuberance'][row+1]], [Cd_low, Cd_up])
            
    # print(Cd)
            
    
    return Cd

def thrust_interp(time, thrust_csv):
    # define thrust for given times
    thrust = 0
        
    for row in range(len(thrust_csv['Time (s)'])-1):

        if thrust_csv['Time (s)'][row] <= time <= thrust_csv['Time (s)'][row+1]:
            thrust = np.interp(time, [thrust_csv['Time (s)'][row], thrust_csv['Time (s)'][row+1]], [thrust_csv['Thrust (N)'][row], thrust_csv['Thrust (N)'][row+1]])

    return thrust