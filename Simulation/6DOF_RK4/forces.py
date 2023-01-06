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
import util.vectors as vct

# Define Objects

class Forces:
    # Define class objects
    atm = atmosphere.Atmosphere()
    motor = motor.Motor() #TODO: Add paramters of motor being used

    def __init__(self):
        pass

    def get_force(self, x_state, flap_ext, tStamp) -> np.ndarray:
        # TODO: Add random disturbances
        # print("State: ", x_state)
        return self.motor.get_thrust(tStamp) + self.aerodynamic_drag_force(x_state, flap_ext)\
        + self.gravitational_force(x_state[0,0], tStamp) + self.wind_force(x_state[0,0], tStamp)

    def get_Cd(self, flap_ext) -> float:
        # TODO: account for area change of flaps in C_d calculation
        return prop.C_d
    
    def aerodynamic_drag_force(self, x_state, flap_ext) -> np.ndarray:
        print("State_A: ", x_state)
        z = x_state[0,0]
        vel = x_state[1].copy()
        v_mag = np.linalg.norm(vel)
        density = self.atm.get_density(z)
        return -vct.norm(vel) * 0.5*density*(v_mag**2)*self.get_Cd(flap_ext)*(prop.A)
    
    def gravitational_force(self, altitude, tStamp) -> np.ndarray:
        # F = GMm/r^2
        total_mass = prop.rocket_dry_mass + self.motor.get_mass(tStamp) #Adding dry mass + motor mass
        return -np.array([(prop.G*prop.m_e*total_mass)/((prop.r_e+altitude)**2), 0, 0])
    
    def wind_force(self, altitude, tStamp) -> np.ndarray:
        wind_vector = self.atm.get_wind_vector(tStamp)
        density = self.atm.get_density(altitude)
        # wind_vector_mag = wind_vector.norm()
        wind_vector_mag = np.linalg.norm(wind_vector)
        # TODO: Change C_d to C_d from the side for wind
        wind_norm = np.zeros((1,3))
        if (wind_vector_mag != 0):
            wind_norm = wind_vector/wind_vector_mag
        return wind_norm * 0.5*density*(wind_vector_mag**2)*prop.C_d*(prop.A_s)

