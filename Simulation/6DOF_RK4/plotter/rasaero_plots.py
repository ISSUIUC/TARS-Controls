import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import math

data = pd.read_csv(os.path.join(os.path.dirname(__file__), '../lookup', 'RASAero.csv'), usecols=['Mach Number', 'Protuberance (%)', 'Alpha (deg)', 'CA Power-Off'])

data = data[data['Alpha (deg)'] == 0.0]
data = data[data['Protuberance (%)'] == 0.0]
print(data.head())

# Generate polynomial fit
c_polyfit = np.polyfit(data['Mach Number'], data['CA Power-Off'], 10)
c_polyfit = np.flip(c_polyfit)
fa_polyfit = np.array([np.sum([c_polyfit[i]*x**i for i in range(len(c_polyfit))]) for x in data['Mach Number']])

def f_true(x):
    x = round(x,2)
    # print(x)

    return data[data['Mach Number'] == x]['CA Power-Off'].values[0]

# Local cubic spline interpolation
def approximate_cubicspline(a, b, n, x_test):
    """Returns function approximations and error for natural cubic spline interpolation.

    This function assumes an f_true(x) function is globally available for calculating the true function value at x.
    
    Parameters
    ----------
    a : float_like
        Lower bound of domain (inclusive)
    b : float_like
        Upper bound of domain (inclusive)
    n : integer
        Number of local splines
    x_test : array_like
        List of inputs to evaluate the approximated function over
        
    Returns
    -------
    fa : array_like
        Vector of approximated function values for each x in x_test
    f : array_like
        Vector of true function
    err : float_like
        Error calculated as a 2-norm using x_test
    x_interpolate : array_like
        List of interpolation points used
    
    """
    
    # get interpolation points (uniform) and delta x
    x_interpolate = np.linspace(a, b, n+1).tolist()
    dx = x_interpolate[1] - x_interpolate[0]

    # get A matrix and f vector
    A = np.zeros((4 * n, 4 * n))
    f = np.zeros(4 * n)
    for i in range(n): # loop through each local spline (i)

        # get first row/column index associated with the ith spline
        ind = i * 4
        
        # update values from condition (1), (5)
        A[ind, ind] = dx**2/6
        A[ind, ind + 1] = 0
        A[ind, ind + 2] = x_interpolate[i]
        A[ind, ind + 3] = 1

        f[ind] = f_true(x_interpolate[i])

        # update values from condition (2), (6)
        A[ind + 1, ind] =  0
        A[ind + 1, ind + 1] = dx**2/6
        A[ind + 1, ind + 2] = x_interpolate[i+1]
        A[ind + 1, ind + 3] = 1

        f[ind + 1] = f_true(x_interpolate[i+1])
        
        if i == n - 1:
            # update values from "extra" condition (7)
            A[ind + 2, 0] = 1

            # update values from "extra" condition (8)
            A[ind + 3, ind + 1] = 1
        else:
            # update values from condition (3)
            A[ind + 2, ind] = 0
            A[ind + 2, ind + 1] = dx/2
            A[ind + 2, ind + 2] = 1
            A[ind + 2, ind + 3] = 0
            A[ind + 2, ind + 4] = dx/2
            A[ind + 2, ind + 5] = 0
            A[ind + 2, ind + 6] = -1
            A[ind + 2, ind + 7] = 0

            # update values from condition (4)
            A[ind + 3, ind] = 0
            A[ind + 3, ind + 1] = 1
            A[ind + 3, ind + 4] = -1

        # update values from conditions (3), (4) and "extra" conditions (7)-(8)
        f[ind + 2] = 0 # f_true(0)
        f[ind + 3] = 0 # f_true(0)
        

    # solve matrix system
    c = np.linalg.solve(A, f)
    # get fa (vector of approximated function values for x_test)
    fa = []
    for x in x_test:
        # get spline index i for x
        if x == x_interpolate[-1]:
            i = n - 1   # index for last interpolation point
        elif x in x_interpolate:
            i = x_interpolate.index(x)  # index for interpolation points
        else:
            i = math.floor(x*(n/3))  # index for test points between interpolation points

        # get first row index associated with the ith spline in c
        ind = i * 4
        
        # get spline i output ("\"" is a line continuation character)
        fa_val = c[ind]/(6*(x_interpolate[i]-x_interpolate[i+1]))*(x-x_interpolate[i+1])**3 + \
                    c[ind+1]/(6*(x_interpolate[i+1]-x_interpolate[i]))*(x-x_interpolate[i])**3 + \
                    c[ind+2]*x + c[ind+3]
        fa.append(fa_val)

    # get f (vector of true function values for x_test)
    f = np.array([f_true(x) for x in x_test])
    # get e (error vector)
    e = f - fa
    # calculate error (2-norm)
    err = np.sqrt(e.T@e)
    return fa, f, err, x_interpolate

def approximate_cubicspline_custom_interpolation(a, b, x_interpolate, x_test):
    """Returns function approximations and error for natural cubic spline interpolation.

    This function assumes an f_true(x) function is globally available for calculating the true function value at x.
    
    Parameters
    ----------
    a : float_like
        Lower bound of domain (inclusive)
    b : float_like
        Upper bound of domain (inclusive)
    n : integer
        Number of local splines
    x_test : array_like
        List of inputs to evaluate the approximated function over
        
    Returns
    -------
    fa : array_like
        Vector of approximated function values for each x in x_test
    f : array_like
        Vector of true function
    err : float_like
        Error calculated as a 2-norm using x_test
    x_interpolate : array_like
        List of interpolation points used
    
    """
    
    # get interpolation points (uniform) and delta x
    # x_interpolate = np.linspace(a, b, n+1).tolist()
    n = len(x_interpolate) - 1

    # get A matrix and f vector
    A = np.zeros((4 * n, 4 * n))
    f = np.zeros(4 * n)
    for i in range(n): # loop through each local spline (i)
        dx = x_interpolate[i + 1] - x_interpolate[i]
        # get first row/column index associated with the ith spline
        ind = i * 4
        
        # update values from condition (1), (5)
        A[ind, ind] = dx**2/6
        A[ind, ind + 1] = 0
        A[ind, ind + 2] = x_interpolate[i]
        A[ind, ind + 3] = 1

        f[ind] = f_true(x_interpolate[i])

        # update values from condition (2), (6)
        A[ind + 1, ind] =  0
        A[ind + 1, ind + 1] = dx**2/6
        A[ind + 1, ind + 2] = x_interpolate[i+1]
        A[ind + 1, ind + 3] = 1

        f[ind + 1] = f_true(x_interpolate[i+1])
        
        if i == n - 1:
            # update values from "extra" condition (7)
            A[ind + 2, 0] = 1

            # update values from "extra" condition (8)
            A[ind + 3, ind + 1] = 1
        else:
            # update values from condition (3)
            A[ind + 2, ind] = 0
            A[ind + 2, ind + 1] = dx/2
            A[ind + 2, ind + 2] = 1
            A[ind + 2, ind + 3] = 0
            A[ind + 2, ind + 4] = dx/2
            A[ind + 2, ind + 5] = 0
            A[ind + 2, ind + 6] = -1
            A[ind + 2, ind + 7] = 0

            # update values from condition (4)
            A[ind + 3, ind] = 0
            A[ind + 3, ind + 1] = 1
            A[ind + 3, ind + 4] = -1

        # update values from conditions (3), (4) and "extra" conditions (7)-(8)
        f[ind + 2] = 0 # f_true(0)
        f[ind + 3] = 0 # f_true(0)
        

    # solve matrix system
    c = np.linalg.solve(A, f)
    # get fa (vector of approximated function values for x_test)
    fa = []
    for x in x_test:
        # get spline index i for x
        if x == x_interpolate[-1]:
            i = n - 1   # index for last interpolation point
        elif x in x_interpolate:
            i = x_interpolate.index(x)  # index for interpolation points
        else:
            i = np.searchsorted(x_interpolate, x) - 1   # index for test points between interpolation points

        # get first row index associated with the ith spline in c
        ind = i * 4
        
        # get spline i output ("\"" is a line continuation character)
        fa_val = c[ind]/(6*(x_interpolate[i]-x_interpolate[i+1]))*(x-x_interpolate[i+1])**3 + \
                    c[ind+1]/(6*(x_interpolate[i+1]-x_interpolate[i]))*(x-x_interpolate[i])**3 + \
                    c[ind+2]*x + c[ind+3]
        fa.append(fa_val)

    # get f (vector of true function values for x_test)
    f = np.array([f_true(x) for x in x_test])
    # get e (error vector)
    e = f - fa
    # calculate error (2-norm)
    err = np.sqrt(e.T@e)
    return fa, f, err, x_interpolate

fa_10, f, err, x_interpolate = approximate_cubicspline(0.01, 3, 10, data['Mach Number'])
fa_15, f, err, x_interpolate = approximate_cubicspline(0.01, 3, 15, data['Mach Number'])
fa_20, f, err, x_interpolate = approximate_cubicspline(0.01, 3, 20, data['Mach Number'])
fa_30, f, err, x_interpolate = approximate_cubicspline(0.01, 3, 30, data['Mach Number'])
fa_300, f, err, x_interpolate = approximate_cubicspline(0.01, 3, 300, data['Mach Number'])
x_interp = [0.01, 0.13, 0.2, 0.4, 0.71, 0.8, 1.06, 1.3, 1.5, 1.7, 2.0, 2.3, 3.0]
fa_custom, f, err, x_interpolate = approximate_cubicspline_custom_interpolation(0.01, 3, x_interp, data['Mach Number'])
plt.plot(data['Mach Number'], data['CA Power-Off'], label='RASAero')
plt.plot(data['Mach Number'], fa_20, label='Cubic Spline, n=20')
plt.plot(data['Mach Number'], fa_30, label='Cubic Spline, n=30')
plt.plot(data['Mach Number'], fa_300, label='Cubic Spline, n=300')
plt.legend()
plt.figure()

plt.plot(data['Mach Number'], (data['CA Power-Off'] - fa_20)/data['CA Power-Off'], label='Error, n=20')
plt.plot(data['Mach Number'], (data['CA Power-Off'] - fa_30)/data['CA Power-Off'], label='Error, n=30')
plt.plot(data['Mach Number'], (data['CA Power-Off'] - fa_300)/data['CA Power-Off'], label='Error, n=300')
plt.legend()
plt.show()
