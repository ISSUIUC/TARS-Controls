'''
Forces on rocket:
- gravity
- motor thrust
- aerodynamic
- wind (comes from atmosphere class)
'''
import motor
import atmosphere
import properties as prop
import numpy as np

# Define Objects

class Forces:
    atm = atmosphere.Atmosphere()

    def __init__(self):
        pass

    def get_force(self, x_state, flap_ext, tStamp):
        return motor.get_thrust + self.aerodynamic_drag_force(x_state, flap_ext) + self.gravitational_force(x_state[0]) + self.wind_force()

    def get_Cd(flap_ext):
        # account for area change of flaps
        return 1
    
    def aerodynamic_drag_force(self, x_state, flap_ext):
        z = x_state[0]
        vel = x_state[1]
        v_mag = np.norm(vel)
        density = self.atm.get_density(z)
        return -(vel/v_mag) * 0.5*density*(v_mag**2)*prop.C_d*(prop.A)
    
    def gravitational_force(altitude):
        # F = GMm/r^2
        return -np.array(0, [(prop.G*prop.m_e*prop.mass)/((prop.r_e+altitude)**2)], 0)
    
    def wind_force(self, altitude, tStamp):
        wind_vector = self.atm.get_wind_vector(tStamp)
        density = self.atm.get_density(altitude)
        wind_vector_mag = wind_vector.norm()
        # Change C_d to C_d from the side for wind
        return (wind_vector/wind_vector_mag) * 0.5*density*(wind_vector_mag**2)*prop.C_d*(prop.A_s)


