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
        
    # convert the angles from degrees to radians
    theta = np.radians(theta - 90)
    phi = np.radians(phi)
    
    # Radius of the earth (km)
    a = 6371.2
    
    # Schmidt quasi-normalized coefficient (are good until January 1, 2025)
    IGRF_lookup = pd.read_csv('Simulation/6DOF_RK4/LookUp/IGRF13coeffs.csv')
    IGRF_dict = IGRF_lookup.set_index(['n','m', 'g/h']).T.to_dict('list')
    #print(IGRF_dict)
    
    '''
    Define the max m and n value (13 currently)
    The max value may change depending on which year of the IGRF model is being used
    To determine the max value, print IGRF_dif and see the very last entry for max values
    '''
    max_M_N = 13
    g = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of g values accessed by [n,m]
    h = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of h values accessed by [n,m]
    
    # Populate g and h
    for key in IGRF_dict:
        if (key[2] == 'g'):
            g[key[0], key[1]] = IGRF_dict.get(key)[0] + (IGRF_dict.get(key)[1] * days / 365)
        else:
            h[key[0], key[1]] = IGRF_dict.get(key)[0] + (IGRF_dict.get(key)[1] * days / 365)
    
    # Initializing Magnetic field values
    Br, Bt, Bp  = 0, 0, 0 # Values in radial direction, theta direction, and phi direction respectively
    
    # Initializing Legendre polynomial values
    # P10 means P(n - 1, m - 0) for naming purposes (same applies for dP)
    P11 = P10 = 1 # Legendre polymial value at cos(theta)
    dP11 = dP10 = 0 # Partial derivative w respect to theta
    
    # figuring this out rn
    for m in range(0, max_M_N + 1):
        for n in range(1, max_M_N + 1):
            if m <= n:
                if n == 1:
                    P2 = np.cos(theta)*P10
                    dP2 = dP10 * np.cos(theta) - P10 * np.sin(theta)*P10
                    P20, dP20 = P10, dP10
                    P10, dP10 = P2, dP2
                elif n == m:
                    P2 = P11 * np.sin(theta)
                    dP2 = dP11 * np.sin(theta) + P11 * np.cos(theta)
                    P11, dP11 = P2, dP2
                    P10, dP10 = P11, dP11
                    P20, dP20 = 0, 0
                else:
                    K = ((n - 1)**2 - m**2) / ((2 * n - 1) * (2 * n - 3))
                    P2 = P10 * np.cos(theta) - K * P20
                    dP2 = dP10 * np.cos(theta) - P10 * np.sin(theta) - K * dP20
                    P20, dP20 = P10, dP10
                    P10, dP10 = P2, dP2
                
                # Calculate Magnetic field values
                Br = Br + (n+1)*((a/r)**(n + 2)) * (g[n, m]*np.cos(m*phi) + h[n,m]*np.sin(m*phi)*P2)
                Bt = Bt + ((a/r)**(n + 2)) * (g[n, m]*np.cos(m*phi) + h[n,m]*np.sin(m*phi)*dP2)
                Bp = Bp + ((a/r)**(n + 2)) * (m*(-g[n, m]*np.cos(m*phi) + h[n,m]*np.sin(m*phi)*P2))
    return Br, -Bt, -Bp / np.sin(theta)

print(magnet(500000, 33, 130, 200))