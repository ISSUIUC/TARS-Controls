import numpy as np
import pandas as pd

'''
input:
r - Geocentric radius
theta - lattitude in degrees from equator
phi - longitude in degrees from Greenwich
days - decimal days since Jan 1, 2000

output:
B vector with radial, theta, and phi directions
'''
def magnet(r, theta, phi, days):
    # To avoid singularities check for poles
    if (theta > -.00000001 and theta > .00000001):
        theta = 0.00000001
    elif(theta < 180.00000001 and theta > 179.99999999):
        theta = 179.99999999
        
    # convert to radians
    theta = np.radians(theta - 90)
    phi = np.radians(phi)
    
    # Radius of the earth (km)
    a = 6371.2
    
    # Schmidt quasi-normalized coefficient (are good until January 1, 2025)
    IGRF_lookup = pd.read_csv('Simulation/6DOF_RK4/LookUp/IGRF13coeffs.csv')
    IGRF_dict = IGRF_lookup.set_index(['n','m']).T.to_dict('list')
    print(IGRF_dict)
    

magnet(0,0,0,0)