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
    motor = motor.Motor()

    def __init__(self):
        pass

    def get_force(self, x_state, flap_ext, time_stamp) -> np.ndarray:
        '''
        Calculates net force felt by rocket while accounting for thrust, drag, gravity, wind

        Args:
            x_state (np.array): State Vector (3x6)
            flap_ext (float): current flap extention config
            time_stamp (float): current time stamp of rocket in simulation
        
        Returns:
            (np.array): 2D array of forces and moments --> ([Fx, Fy, Fz], [Mx, My, Mz])
        '''
        # TODO: Add random disturbances
        thrust = self.motor.get_thrust(time_stamp)
        drag = self.aerodynamic_drag_force(x_state, flap_ext)
        grav = self.gravitational_force(x_state[0,0], time_stamp)
        wind = self.wind_force(x_state[0,0], time_stamp)
        # print(self.motor.get_thrust(time_stamp))
        force = thrust + drag + wind + grav
        moment = np.cross(-prop.cm, force)
        return np.array([force, moment])

    def get_Cd(self, flap_ext) -> float:
        # TODO: account for area change of flaps in C_d calculation
        '''
        References lookup table to find C_d based on flap extension

        Args:
            flap_ext (float): current flap extention config
        
        Returns:
            (float): coefficient of drag based on current config of flaps
        '''
        return prop.C_d
    
    def aerodynamic_drag_force(self, x_state, flap_ext) -> np.ndarray:
        '''
        Calculates aerodynamic drag force acting on rocket based on velocity and altitude

        Args:
            x_state (np.array): State Vector (3x6)
            flap_ext (float): current flap extention config
        
        Returns:
            (np.array): vector of aerodynamic forces in each axis (1x3)
        '''
        z = x_state.copy()[0,0]
        vel = x_state[1].copy()
        v_mag = np.linalg.norm(vel)
        density = self.atm.get_density(z)
        return -vct.norm(vel) * 0.5*density*(v_mag**2)*self.get_Cd(flap_ext)*(prop.A)
    
    def gravitational_force(self, altitude, time_stamp) -> np.ndarray:
        '''
        Calculates gravitational force acting on rocket based on altitude
        Relevant Equations:
            F = GMm/r^2

        Args:
            altitude (float): current altitude of rocket
            time_stamp (float): current time stamp of rocket in simulation
        
        Returns:
            (np.array): vector of gravitational forces on each axis (1x3)
        '''
        total_mass = prop.rocket_dry_mass + self.motor.get_mass(time_stamp) #Adding dry mass + motor mass
        return -np.array([(prop.G*prop.m_e*total_mass)/((prop.r_e+altitude)**2), 0, 0])
    
    def wind_force(self, altitude, time_stamp) -> np.ndarray:
        '''
        Calculates wind force acting on rocket based on determined wind and altitude

        Args:
            altitude (float): current altitude of rocket
            time_stamp (float): current time stamp of rocket in simulation
        
        Returns:
            (np.array): vector of wind forces on each axis (1x3)
        '''
        wind_vector = self.atm.get_wind_vector(time_stamp)
        density = self.atm.get_density(altitude)
        wind_vector_mag = np.linalg.norm(wind_vector)
        # TODO: Change C_d to C_d from the side for wind
        wind_norm = vct.norm(wind_vector)
        return wind_norm * 0.5*density*(wind_vector_mag**2)*prop.C_d*(prop.A_s)

