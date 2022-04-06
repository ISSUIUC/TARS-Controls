import numpy as np
import pandas as pd
import csv
import sys
from scipy.integrate import quad

def find_prop_mass(time):

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