'''
References:
https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html - IGRF Coefficients
https://hanspeterschaub.info/Papers/UnderGradStudents/MagneticField.pdf - MATLAB implementation used as a guide
https://geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html - Online Magnetic field calculator for testing
'''


import numpy as np
import pandas as pd

def magnet(r, theta, phi, days):
    '''Magentic field model for the Earth's magnetic field

    Args:
        r: Geocentric radius
        theta: lattitude in degrees (North) from equator
        phi: longitude in degrees (East) from Greenwich
        days: decimal days since Jan 1, 2020

    Return:
        B vector in x, y, and z directions respectively in nT
    '''
    # To avoid singularities check for poles
    if (theta > -.00000001 and theta < .00000001):
        theta = 0.00000001
    elif(theta < 180.00000001 and theta > 179.99999999):
        theta = 179.99999999
        
    # convert the angles from degrees to radians
    theta = np.radians(90 - theta)
    phi = np.radians(phi)
    
    # Radius of the earth (km)
    a = 6371.2
    
    # Schmidt quasi-normalized coefficient (are good until January 1, 2025)
    gS, hS, max_M_N = quasi_normalize('Simulation/6DOF_RK4/LookUp/IGRF13coeffs.csv')
    
    
    g = np.zeros((max_M_N + 1, max_M_N + 1))
    h = np.zeros((max_M_N + 1, max_M_N + 1))
    for i in range(len(gS[:])):
        g[int(gS[i, 0]), int(gS[i, 1])] = gS[i, 2] + gS[i, 3]*days/365
        h[int(hS[i, 0]), int(hS[i, 1])] = hS[i, 2] + hS[i, 3]*days/365
    
    # Initializing Magnetic field values
    Br, Bt, Bp  = 0, 0, 0 # Values in radial direction, theta direction, and phi direction respectively
    
    # Initializing Legendre polynomial values
    # P10 means P(n - 1, m - 0) for naming purposes (same applies for dP)
    P11 = P10 = 1 # Legendre polymial value at cos(theta)
    dP11 = dP10 = 0 # Partial derivative w respect to theta
    
    # Recursively finding legendre polynomial and magnetic field vals
    for m in range(max_M_N + 1):
        for n in range(1, max_M_N + 1):
            if m <= n:
                # Calculate Legendre polynomials recursively
                # P_curr and dP_curr corresponds to current n and m value
                if n == m:
                    P_curr = P11 * np.sin(theta)
                    dP_curr = dP11 * np.sin(theta) + P11 * np.cos(theta)
                    P11, dP11 = P_curr, dP_curr
                    P10, dP10 = P11, dP11
                    P20, dP20 = 0, 0
                    
                elif n == 1:
                    P_curr = np.cos(theta)*P10
                    dP_curr = dP10 * np.cos(theta) - np.sin(theta)*P10
                    P20, dP20 = P10, dP10
                    P10, dP10 = P_curr, dP_curr
                    
                else:
                    K = ((n - 1)**2 - m**2) / ((2*n - 1) * (2*n - 3))
                    P_curr = P10 * np.cos(theta) - K * P20
                    dP_curr = dP10 * np.cos(theta) - P10 * np.sin(theta) - K * dP20
                    P20, dP20 = P10, dP10
                    P10, dP10 = P_curr, dP_curr
                
                # Calculate Magnetic field values
                Br = Br + ((a/r)**(n + 2)) * (n+1) * ((g[n,m]*np.cos(m*phi) + h[n, m]*np.sin(m*phi)) * P_curr)
                Bt = Bt + ((a/r)**(n + 2)) * ((g[n, m]*np.cos(m*phi) + h[n,m]*np.sin(m*phi))*dP_curr)
                Bp = Bp + ((a/r)**(n + 2)) * (m*(-g[n, m]*np.sin(m*phi) + h[n,m]*np.cos(m*phi))*P_curr)
    
    return sphere_to_cartesian(Br, -Bt, -Bp / np.sin(theta))

def sphere_to_cartesian(Br, Bt, Bp):
    '''Used to convert spherical magnetic field to cartesian coordinates
    
    Args:
        Br (float): B vector radial direction
        Bt (float): B vector theta direction
        Bp (float): B vector phi direction
        
    Returns: (all directions by convention)
        Bx (float): Magnetic field in North direction
        By (float): Magnetic field in East Direction
        Bz (float): Magnetic field in Downward Direction
    '''
    e = 0
    Bx = -Bt*np.cos(e) - Br*np.sin(e)
    By = Bp
    Bz = Bt*np.sin(e) - Br*np.cos(e)
    return Bx, By, Bz


def quasi_normalize(file_lookup_path):
    '''Quasi-normalizing the input coefficient file

    Args:
        file_lookup_path: Path to the IGRF coefficient file
        
    Returns: 
        gs: Quasi-normalized g values
        hs: Quasi-normalized h values
        max_M_N: Max value of m and n
    '''
    IGRF_lookup = pd.read_csv(file_lookup_path)
    IGRF_dict = IGRF_lookup.set_index(['n','m', 'g/h']).T.to_dict('list')
    
    '''
    Define the max m and n value (13 currently)
    The max value may change depending on which year of the IGRF model is being used
    To determine the max value, print IGRF_dict and see the very last entry for max values
    There's probably a better way to do this, but I didn't want to add another loop in here
    '''
    max_M_N = 13
    G = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of g values accessed by [n,m]
    H = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of h values accessed by [n,m]
    GSV = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of G coeff SV values accessed by [n,m]
    HSV = np.zeros((max_M_N + 1, max_M_N + 1)) # 2d array of H coeff SV values accessed by [n,m]
    
    # Populate G, H, HSV, HSV with corresponding non-normalized coeffs and SV's
    for key in IGRF_dict:
        if (key[2] == 'g'):
            G[key[0], key[1]] = IGRF_dict.get(key)[0]
            GSV[key[0], key[1]] = IGRF_dict.get(key)[1]
        else:
            H[key[0], key[1]] = IGRF_dict.get(key)[0]
            HSV[key[0], key[1]] = IGRF_dict.get(key)[1]
    
    # Finding the size gS and hS will be
    size_int = 0
    for n in range(max_M_N + 1):
        for m in range(n + 1):
            size_int += 1
    
    # Initializing S, gS, and hS
    S = np.zeros((max_M_N + 1, max_M_N + 1))
    gS = np.zeros((size_int, 4))
    hS = np.zeros((size_int, 4))
    
    count = 0
    for n in range(1, max_M_N + 1):
        for m in range(n + 1):
            if m > 1:
                S[n, m] = S[n, m - 1]*((n - m + 1)/(n + m)) ** 0.5
            elif m > 0:
                S[n, m] = S[n, m - 1]*(2 * (n - m + 1)/(n + m)) ** 0.5
            elif n == 1:
                S[n, 0] = 1
            else:
                S[n, 0] = S[n - 1, 0]*(2*n-1)/n
            gS[count] = [n, m, G[n, m] * S[n, m], GSV[n, m] * S[n, m]]
            hS[count] = [n, m, H[n, m] * S[n, m], HSV[n, m] * S[n, m]]
            
            count += 1
    return gS, hS, max_M_N
