from numpy import pi


def ft_to_m(measurement): 
    return (measurement / 3.2808) 

def c_to_k(measurement):
    return (measurement + 273.15)

def deg_to_rad(measurement):
    return (measurement*pi/180)

def rad_to_deg(measurement):
    return (measurement*180/pi)