'''
For each unit test, the code calculations as of 9/9/23 is run by the inner test function. 
This is saved to the outer test function which is not run by the inner test function. If
the test output by the test method is identical to the test output by a modified function 
in the main dynamics code, true will be returned and the test will have passed. Else, false
is returned.

TODO: 
- create tests for all "units"
- enter test values
- create naming conventions to avoid double-naming
- create function running every test and returning false if any test fails and if the test fails, which tests failed
    -also maybe why??
- verify current code calcs
'''


import dynamics.motor as motor
import environment.atmosphere as atmosphere
import properties.properties as prop
import numpy as np
import util.vectors as vct
import pandas as pd
import os

class ForceTest():

    def __init__(self, max_ext_length, cm, cp, A, A_s, rocket_dry_mass, motor, atm): # Do we need this?
        self.max_ext_length = max_ext_length
        self.cm = cm
        self.cp = cp
        self.A = A
        self.A_s = A_s
        self.rocket_dry_mass = rocket_dry_mass
        self.motor = motor
        self.atm = atm

    def test_forces(self, x_state, flap_ext, time_stamp, density_noise=False, test) -> np.ndarray:
        '''Calculates net force felt by rocket while accounting for thrust, drag, gravity, wind

        Args:
            x_state (np.array): State Vector [4x3]
            flap_ext (float): current flap extention config
            time_stamp (float): current time stamp of rocket in simulation
            test ALWAYS is false
        
        Returns:
            Boolean for if test passes or not
        '''
        if(test):  
            # print("State: ", x_state)
            alt = x_state.copy()[0,0]
            density = self.atm.get_density(alt, noise=density_noise, position=x_state[0])
            thrust = self.motor.get_thrust(time_stamp)
            # wind_vector = self.atm.get_nominal_wind_direction() * self.atm.get_nominal_wind_magnitude()
            wind_vector = self.atm.get_wind_vector(time_stamp)
            alpha = self.get_alpha(x_state, wind_vector)
            drag = self.aerodynamic_force(x_state, density, wind_vector, alpha, self.rasaero, thrust.dot(thrust) > 0, flap_ext)
            grav = self.gravitational_force(alt, time_stamp)
            force = vct.body_to_world(*x_state[2],thrust + drag) + grav
            moment = vct.body_to_world(*x_state[2], np.cross(-self.cm, thrust) + self.aerodynamic_moment(drag))
            # print(self.aerodynamic_moment(drag))
            return np.array([force, moment]), alpha
        if not(test):
            forcesTest = test_forces(self, x_state, flap_ext, time_stamp, density_noise=False, True) #put in our test values here
            if(forcesTest == getForce(self, x_state, flap_ext, time_stamp, density_noise=False)):
                return True
            else:
                return False

        