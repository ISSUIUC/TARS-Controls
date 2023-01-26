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
import pandas as pd
import os

# Define Objects

class Forces:
    # Define class objects
    atm = atmosphere.Atmosphere()
    motor = motor.Motor()

    rasaero_file_location = os.path.join(os.path.dirname(__file__), prop.rasaero_lookup_file)
    rasaero = pd.read_csv(rasaero_file_location)

    def __init__(self):
        pass

    def get_force(self, x_state, flap_ext, time_stamp) -> np.ndarray:
        '''
        Calculates net force felt by rocket while accounting for thrust, drag, gravity, wind

        Args:
            x_state (np.array): State Vector [4x3]
            flap_ext (float): current flap extention config
            time_stamp (float): current time stamp of rocket in simulation
        
        Returns:
            (np.array): 2D array of forces and moments --> ([Fx, Fy, Fz], [Mx, My, Mz])
        '''
        # TODO: Add random disturbances
        # print("State: ", x_state)
        alt = x_state.copy()[0,0]
        density = self.atm.get_density(alt)
        thrust = self.motor.get_thrust(time_stamp)
        wind_vector = self.atm.get_nominal_wind_direction() * self.atm.get_nominal_wind_magnitude()
        drag = self.aerodynamic_force(x_state, density, wind_vector, self.rasaero, thrust.dot(thrust) > 0, flap_ext)
        grav = self.gravitational_force(alt, time_stamp)
        force = vct.body_to_world(*x_state[2],thrust + drag) + grav
        moment = vct.body_to_world(*x_state[2], np.cross(-prop.cm, thrust) + self.aerodynamic_moment(drag))
        return np.array([force, moment])

    def get_Ca(self, x_state, wind_vector, rasaero, before_burnout, flap_ext) -> float:
        # TODO: account for area change of flaps in C_a calculation
        '''
        References lookup table to find C_a based on flap extension

        Args:
            flap_ext (float): current flap extention config
        
        Returns:
            (float): coefficient of drag based on current config of flaps
        '''
        alt = x_state[0,0]
        vel = x_state[1]
        mach_number = np.linalg.norm(vel) / self.atm.get_speed_of_sound(alt)
        
        incident_velocity = vct.world_to_body(*x_state[2], vct.norm(vel + wind_vector))
        orientation = vct.world_to_body(*x_state[2], np.array([1,0,0]))
        alpha = np.arccos(np.dot(incident_velocity, orientation))
        # Define mach number for csv lookup, rounded to hundreds place
        mach = round(mach_number, 2)
        
        # round AoA to closest integer
        if(np.linalg.norm(incident_velocity) == 0):
            alpha = 0
        ang_of_att = round(np.degrees(alpha))

        #Define blank upper and lower Cds
        Ca_low = 0
        Ca_up = 0

        # define csv file to search through
        #csv_file = csv.reader(open('RASAero.csv', 'r'))
        csv_file = rasaero

        # Define protuberance percentage of full extension given current extension
        protub_perc = flap_ext/prop.max_ext_length

        # Define starting Ca
        Ca = 0
        
        ca_vals = csv_file["CA Power-Off"]
        if (before_burnout):
            ca_vals = csv_file["CA Power-On"]
        
        # Find indices where the mach values match up
        mach_indices = np.where(csv_file['Mach Number'] == mach)[0]    
        if len(mach_indices) == 0:
            return Ca
            
        # Interpolate to find Cn value
        for idx in range(mach_indices[0], mach_indices[-1] + 1):
            if (ang_of_att == csv_file['Alpha (deg)'][idx] and csv_file['Protuberance (%)'][idx] <= protub_perc <= csv_file['Protuberance (%)'][idx+1]):
                
                Ca_low = ca_vals[idx]
                Ca_up = ca_vals[idx+1]

                Ca = np.interp(protub_perc, [csv_file['Protuberance (%)'][idx], csv_file['Protuberance (%)'][idx+1]], [Ca_low, Ca_up])    
        return Ca
    
    def get_Cn(self, x_state, wind_vector, rasaero, flap_ext) -> float:
        # TODO: account for area change of flaps in C_n calculation
        '''
        References lookup table to find C_n based on flap extension

        Args:
            flap_ext (float): current flap extention config
        
        Returns:
            (float): coefficient of drag based on current config of flaps
        '''
        alt = x_state[0,0]
        vel = x_state[1]
        mach_number = np.linalg.norm(vel) / self.atm.get_speed_of_sound(alt)
        
        incident_velocity = vct.world_to_body(*x_state[2], vct.norm(vel + wind_vector))
        orientation = vct.world_to_body(*x_state[2], np.array([1,0,0]))
        alpha = np.arccos(np.dot(incident_velocity, orientation))
        # Define mach number for csv lookup, rounded to hundreds place
        mach = round(mach_number, 2)
        
        # round AoA to closest integer
        if(np.linalg.norm(incident_velocity) == 0):
            alpha = 0
        
        ang_of_att = round(np.degrees(alpha))
        #Define blank upper and lower Cds
        Cn_low = 0
        Cn_up = 0

        # define csv file to search through
        #csv_file = csv.reader(open('RASAero.csv', 'r'))
        csv_file = rasaero

        # Define protuberance percentage of full extension given current extension
        protub_perc = flap_ext/prop.max_ext_length

        # Define starting Cn
        Cn = 0
        
        cn_vals = csv_file["CN Total"]
        
        # Find indices where the mach values match up
        mach_indices = np.where(csv_file['Mach Number'] == mach)[0]    
        if len(mach_indices) == 0:
            return Cn
            
        # Interpolate to find Cn value
        for idx in range(mach_indices[0], mach_indices[-1] + 1):
            if (ang_of_att == csv_file['Alpha (deg)'][idx] and csv_file['Protuberance (%)'][idx] <= protub_perc <= csv_file['Protuberance (%)'][idx+1]):
                
                Cn_low = cn_vals[idx]
                Cn_up = cn_vals[idx+1]

                Cn = np.interp(protub_perc, [csv_file['Protuberance (%)'][idx], csv_file['Protuberance (%)'][idx+1]], [Cn_low, Cn_up])    
        return Cn

    def aerodynamic_force(self, x_state, density, wind_vector, rasaero, before_burnout, flap_ext) -> np.ndarray:
        '''
        Calculates aerodynamic drag force acting on rocket based on velocity and altitude

        Args:
            x_state (np.array): State Vector [4x3]
            flap_ext (float): current flap extention config
        
        Returns:
            (np.array): vector of aerodynamic forces in each axis [1x3]
        '''
        vel = vct.world_to_body(*x_state[2].copy(), x_state[1].copy() - wind_vector.copy())
        C_a = self.get_Ca(x_state, wind_vector, rasaero, before_burnout, flap_ext)
        C_n = self.get_Cn(x_state, wind_vector, rasaero, flap_ext)

        roll_aero = np.arctan2(x_state[1,2], x_state[1,1])

        C_n_y = np.abs(C_n * np.cos(roll_aero)) #TODO: Check with other values
        C_n_z = np.abs(C_n * np.sin(roll_aero)) #TODO: Check with other values
        aero_force = -0.5*np.array([np.sign(vel[0])*vel[0]**2 * C_a*density*prop.A, 
                                    np.sign(vel[1])*vel[1]**2 * C_n_y*density*prop.A_s, 
                                    np.sign(vel[2])*vel[2]**2 * C_n_z*density*prop.A_s])
        return aero_force
    
    def gravitational_force(self, altitude, time_stamp) -> np.ndarray:
        '''
        Calculates gravitational force acting on rocket based on altitude
        Relevant Equations:
            F = GMm/r^2

        Args:
            altitude (float): current altitude of rocket
            time_stamp (float): current time stamp of rocket in simulation
        
        Returns:
            (np.array): vector of gravitational forces on each axis [1x3]
        '''
        total_mass = prop.rocket_dry_mass + self.motor.get_mass(time_stamp) #Adding dry mass + motor mass
        # return np.array([-9.81*total_mass, 0, 0])
        return -np.array([(prop.G*prop.m_e*total_mass)/((prop.r_e+altitude)**2), 0, 0])

    def aerodynamic_moment(self, aerodynamic_force):
        aerodynamic_moment = np.cross(prop.cp - prop.cm, aerodynamic_force)
        return aerodynamic_moment
