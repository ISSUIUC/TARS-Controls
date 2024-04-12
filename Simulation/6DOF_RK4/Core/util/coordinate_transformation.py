import numpy as np


def sphr2cart(Br, Bt, Bp):
    '''
    Converts magnetic field from spherical coordinates to cartesian coordinates

    Variables:
        e (float): episilon (term used to account for the oblateness of the Earth)

    Args:
        Br (float): Magnetic Field in radial direction
        Bt (float): Magnetic Field in theta direction
        Bp (float): Magnetic Field in phi direction

    Returns:
        Bx (float): Magnetic Field in "North" direction
        By (float): Magnetic Field in "East" direction
        Bz (float): Magnetic Field in down direction
    '''
    e = 0*np.pi/180
    Bx = -Bt*np.cos(e) + Br*np.sin(e)
    By = Bp
    Bz = Bt*np.cos(e) - Br*np.sin(e)
    return (Bx, By, Bz)


def sphr2geocart(Br, Bt, Bp, lst, lat):
    '''
    Converts magnetic field from spherical coordinates to geocentric inertial cartesian coordinates

    Args:
        Br (float): Magnetic Field in radial direction
        Bt (float): Magnetic Field in theta direction
        Bp (float): Magnetic Field in phi direction
        LST (float): Local Sidereel Time (in degrees)
        lat (float): Latitude (in degrees) (positive is north)

    Returns:
        Bx (float): Magnetic Field in "North" direction
        By (float): Magnetic Field in "East" direction
        Bz (float): Magnetic Field in down direction
    '''

    # Convert Degrees to Radians
    lst = lst * np.pi / 180
    lat = lat * np.pi / 180

    Bx = (Br*np.cos(lat)+Bt*np.sin(lat))*np.cos(lst) - Bp*np.sin(lst)
    By = (Br*np.cos(lat)+Bt*np.sin(lat))*np.sin(lst) + Bp*np.cos(lst)
    Bz = (Br*np.sin(lat)+Bt*np.cos(lat))
    return (Bx, By, Bz)
