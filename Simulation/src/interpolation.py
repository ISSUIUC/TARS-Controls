import numpy as np
import pandas as pd
import csv
import src.atmosphere as atmosphere
import sys

# for inputs, extension is a value in inches, between 0 and 0.5. rasaero is the csv file pulled from Lookup

def cd_interpolation(altitude, velocity, ang_of_att, extension, rasaero):
    
    # Use atmosphere in src to find speed of sound given altitude
    a = atmosphere.speed_sound(altitude)

    # Define mach number for csv lookup, rounded to hundreds place
    mach = round((velocity/a), 2)
    print(mach)

    # round AoA to closest integer
    ang_of_att = round(ang_of_att)

    # Define numbers to lookup in csv
    protub_off = 0
    protub_on = 1

    #Define blank upper and lower Cds
    Cd_low = 0
    Cd_up = 0

    # define csv file to search through
    #csv_file = csv.reader(open('RASAero.csv', 'r'))
    csv_file = rasaero

    # for loop to find upper Cd value of interpolation
    for row in range(len(csv_file['mach'])):
        # see if row meets values for the lower Cd

        if (mach == csv_file['mach'][row]) & (ang_of_att == csv_file['alpha_deg'][row]) & (protub_off == csv_file['protuberance'][row]):
            Cd_low = csv_file['cd_power_off'][row]

    print(Cd_low)

    # for loop to find lower Cd value of interpolation
    for row in range(len(csv_file['mach'])):
        # see if row meets values for the lower Cd
        if (mach == csv_file['mach'][row]) & (ang_of_att == csv_file['alpha_deg'][row]) & (protub_on == csv_file['protuberance'][row]):
            Cd_up = csv_file['cd_power_off'][row]

    print(Cd_up)
    # use numpy to interpoalate the Cd using upper and lower bounds pulled from csv and full/no extension
    Cd = np.interp(extension, [0, 0.5], [Cd_low, Cd_up])
    return Cd



def thrust_interp(time, thrust_csv):
    # define thrust for given times
    thrust = 0
    
    for row in range(len(thrust_csv['Time (s)'])-1):

        if thrust_csv['Time (s)'][row] <= time <= thrust_csv['Time (s)'][row+1]:
            thrust = np.interp(time, [thrust_csv['Time (s)'][row], thrust_csv['Time (s)'][row+1]], [thrust_csv['Thrust (N)'][row], thrust_csv['Thrust (N)'][row+1]])

    return thrust