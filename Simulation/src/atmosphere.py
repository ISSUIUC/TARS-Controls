# ACCURATE UP TO 11000m
def pressure(z):
    # temperature under standard condition (15 degrees C at sealevel) kelvin
    P_0 = 101325

    # pressure under standard condition in (Pa)
    T_0 = 288.16 

    # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
    b = 0.0065

    # gravitational constant 
    g = 9.81

    R = 287.05

    return P_0 * ((T_0 +(z)*b)/T_0)**(-g/(b*R))

def alt_from_pressure(p):
    # temperature under standard condition (15 degrees C at sealevel) kelvin
    P_0 = 101325

    # pressure under standard condition in (Pa)
    T_0 = 288.16 

    # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
    b = 0.0065

    # gravitational constant 
    g = 9.81

    R = 287.05

    pressureRatio = p/P_0
    return -(T_0*((pressureRatio)**(b*R/(g)) - 1) * (pressureRatio)**(-b*R/(g)))/b

def density(z):
    
    # temperature under standard condition (15 degrees C at sealevel) kelvin
    T_0 = 288.16 

    # pressure under standard condition in (Pa)
    P_0 = 101325

    # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
    b = 0.0065
    
    # gravitational constant 
    g = 9.81

    # air density under standard condition (kg/m3)
    rho_0 = 1.225 
    
    # ideal gas constant (J/kg K)
    R = 287.05

    # Return the density
    return rho_0*(1 - ((b*z)/T_0))**(g/(R*b))*(T_0/(T_0 - b*z))

def viscosity(z):
    
    #TODO: What are these constants?
    u0 = 1.716e-5 
    Su = 111
    
    # Reference temperature for Sutherland's Law
    T0 = 273.15
    
    # temperature under standard condition (15 degrees C at sealevel) kelvin
    T_0 = 288.16
    
    # Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
    b = 0.0065
    
    # Current temperature based on altitute 
    T  = T_0 - b*z
    
    # Return the dynamic viscosity
    return ((T/T0)**(1.5))*((T0 + Su)/(T + Su))*u0

def speed_sound(altitude):
    #* the first 11000 m of flight, the temperature varies linearly with altitude (k/m)
    dt_dh = 6.5e-3 
    #* specific heat ratio of air 
    r = 1.4 
    #* ideal gas constant of air R 
    R = 287.05
    #* temperature under standard condition (15 degrees C at sealevel) kelvin
    T_0 = 288.16 
    #* Current temperature based on altitute 
    T = T_0 - dt_dh*altitude
    #* Calculating speed of sound 
    a_H = (r*R*T)**0.5

    return a_H

