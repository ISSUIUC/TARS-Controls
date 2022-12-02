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
    motor = motor.Motor()

    def __init__(self):
        pass

    def get_force(self, x_state, flap_ext, tStamp) -> float:
        return self.motor.get_thrust(tStamp) + self.aerodynamic_drag_force(x_state, flap_ext) + self.gravitational_force(x_state[0], tStamp) + self.wind_force()

    def get_Cd(flap_ext) -> float:
        # TODO: account for area change of flaps in C_d calculation
        return prop.C_d
    
    def aerodynamic_drag_force(self, x_state, flap_ext) -> float:
        z = x_state[0]
        vel = x_state[1]
        v_mag = np.norm(vel)
        density = self.atm.get_density(z)
        return -(vel/v_mag) * 0.5*density*(v_mag**2)*self.get_Cd(flap_ext)*(prop.A)
    
    def gravitational_force(self, altitude, tStamp) -> float:
        # F = GMm/r^2
        total_mass = prop.mass + self.motor.get_mass(tStamp) #Adding dry mass + motor mass
        return -np.array(0, [(prop.G*prop.m_e*total_mass)/((prop.r_e+altitude)**2)], 0)
    
    def wind_force(self, altitude, tStamp) -> float:
        wind_vector = self.atm.get_wind_vector(tStamp)
        density = self.atm.get_density(altitude)
        wind_vector_mag = wind_vector.norm()
        # TODO: Change C_d to C_d from the side for wind
        return (wind_vector/wind_vector_mag) * 0.5*density*(wind_vector_mag**2)*prop.C_d*(prop.A_s)


