import numpy as np
import pandas as pd
import csv
import sys
from scipy.integrate import quad

# output propellant mass expended at a given time, using the best fit equation for the April launch motor

def find_prop_mass_april(time):

    def m_dot_eqn(t):
        
        # from motor datasheet
        Isp = 209 # s
        g = 9.81 # m/s^2

        # mass flow equation is m_dot = T/(Isp g) where T is thrust
        # below uses a fourth order best-fit for the thrust curve to approxiamnte flow rate
        # Motor is Aerotech N2500T (April Launch)
        return (1/(Isp*g))*(24.878*(t**4) - 322.81*(t**3) + 900.96*(t**2) - 677.87*(t) + 2877.6)


    # m_prop is mass of propellant consumed at given time, err is the expected error from true solution
    m_prop, err = quad(m_dot_eqn, 0.082, time)

    return m_prop


# output propellant mass expended at a given time, using the best fit equation for the IREC launch motor

def find_prop_mass_irec(time):

    def m_dot_eqn(t):
        
        # from motor datasheet
        Isp = 228 # s
        g = 9.81 # m/s^2

        # mass flow equation is m_dot = T/(Isp g) where T is thrust
        # below uses a fourth order best-fit for the thrust curve to approxiamnte flow rate
        # Motor is Aerotech N2500T (April Launch)
        return (1/(Isp*g))*(134.96*(t**6) - 1319.1*(t**5) + 4701.6*(t**4) - 7829*(t**3) + 6038.1*(t**2) - 1532.1*(t) + 6733.6)


    # m_prop is mass of propellant consumed at given time, err is the expected error from true solution
    m_prop, err = quad(m_dot_eqn, 0.019, time)

    return m_prop