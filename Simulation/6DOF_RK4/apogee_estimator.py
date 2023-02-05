import numpy as np
import scipy
import matplotlib.pyplot as plt
import forces
import util.vectors as vct
import properties as prop
import atmosphere
import motor
import os
import pandas as pd

class Apogee: 

    atm = atmosphere.Atmosphere()
    motor = motor.Motor()
    rasaero_file_location = os.path.join(os.path.dirname(__file__), prop.rasaero_lookup_file)
    rasaero = pd.read_csv(rasaero_file_location)


    def __init__(self, state, dt):
        '''
        self.state:
        
        x-axis:
        [[pos],
         [vel],
         [acc]]
        '''
        self.state = state[:3][0].copy()
        self.dt = dt
        self.flap_ext = 0.
    
    def set_params(self, state):
        self.state = state[:3][0].copy()

    def RK4(self):
        k1_v = self.state[2]
        k2_v = self.step_v(self.state[1] + k1_v*(self.dt/2))
        k3_v = self.step_v(self.state[1] + k2_v*(self.dt/2))
        k4_v = self.step_v(self.state[1] + k3_v*(self.dt))

        v = (self.state[1] + (1/6)*(k1_v+(2*k2_v)+(2*k3_v)+k4_v)*self.dt)
    
        k1_p = self.state[1]
        k2_p = self.step_p(self.state[0] + k1_p*(self.dt/2))
        k3_p = self.step_p(self.state[0] + k2_p*(self.dt/2))
        k4_p = self.step_p(self.state[0] + k3_p*(self.dt))

        p = (self.state[0] + (1/6)*(k1_p+(2*k2_p)+(2*k3_p)+k4_p)*self.dt)

        a = self.get_accel()

        self.state = np.array([p,v,a])
    
    def step_p(self, y0, y1, dt):
        return (y1-y0)/dt
    
    def get_Ca(self) -> float:
        # TODO: account for area change of flaps in C_a calculation
        '''
        References lookup table to find C_a based on flap extension

        Args:
            flap_ext (float): current flap extention config
        
        Returns:
            (float): coefficient of drag based on current config of flaps
        '''
        alt = self.state[0]
        vel = self.state[1]
        
        # Define mach number for csv lookup, rounded to hundreds place
        mach_number = np.linalg.norm(vel) / self.atm.get_speed_of_sound(alt)
        mach = round(mach_number, 2)
        
        #Define blank upper and lower Cds
        Ca_low = 0
        Ca_up = 0

        # define csv file to search through
        csv_file = self.rasaero

        # Define protuberance percentage of full extension given current extension
        protub_perc = self.flap_ext/prop.max_ext_length

        # Define starting Ca
        Ca = 0
        
        ca_vals = csv_file["CA Power-Off"]

        # Find indices where the mach values match up
        mach_indices = np.where(csv_file['Mach Number'] == mach)[0]    
        if len(mach_indices) == 0:
            return Ca
            
        # Interpolate to find Cn value
        for idx in range(mach_indices[0], mach_indices[-1] + 1):
            if (0 == csv_file['Alpha (deg)'][idx] and csv_file['Protuberance (%)'][idx] <= protub_perc <= csv_file['Protuberance (%)'][idx+1]):
                
                Ca_low = ca_vals[idx]
                Ca_up = ca_vals[idx+1]

                Ca = np.interp(protub_perc, [csv_file['Protuberance (%)'][idx], csv_file['Protuberance (%)'][idx+1]], [Ca_low, Ca_up])    
        return Ca


    def get_accel(self) -> np.ndarray:
        '''
        Calculates net force felt by rocket while accounting for thrust, drag, gravity, wind
        
        Returns:
            (np.array): 2D array of forces and moments --> ([Fx, Fy, Fz], [Mx, My, Mz])
        '''
        # TODO: Add random disturbances
        # print("State: ", x_state)
        C_a = self.get_Ca()
        alt = self.state[0]
        density = self.atm.get_density(alt)
        drag = -0.5*(self.state[1]**2 * C_a*density*prop.A)
        grav = self.gravitational_force(alt)
        force =  drag + grav
        return force/prop.rocket_dry_mass

    def gravitational_force(self, altitude) -> np.ndarray:
        '''
        Calculates gravitational force acting on rocket based on altitude
        Relevant Equations:
            F = GMm/r^2

        Args:
            altitude (float): current altitude of rocket
        
        Returns:
            (np.array): vector of gravitational forces on each axis [1x3]
        '''
        return -(prop.G*prop.m_e*prop.rocket_dry_mass)/((prop.r_e+altitude)**2)

    def step_v(self):
        
        '''
        Calculates slope of v over given delta t for state propogation
        
        Returns:
            (np.array): rate of change of velocity (acceleration) in form of state vector
        '''
        return self.get_accel()    
    
    def predict_apogee(self, current_state):
        self.set_params(current_state.copy())
        while (self.state[1] > 0):
            self.RK4()

        return self.state[0]


if __name__=='__main__':
    import matplotlib.pyplot as plt
    data = pd.read_csv(os.path.join(os.path.dirname(__file__), "flight_computer_20221029.csv"))
    # print(data.head())
    state_estimate = data["state_est_x"]
    ax = data["ax"]
    rocket_estimated_apogee = data["state_est_apo"]
    timestamps = data["timestamp_ms"]
    # print(state_estimate.head())
    # print(ax.head())
    #TODO: add testing for apogee estimator and fix state input
    apogee_estimator = Apogee(state_estimate[0], 0.01)
    apogee_estimates = np.array([])
    for estimate in state_estimate:
        apogee_estimate = apogee_estimator.predict_apogee(state_estimate[0])
        np.append(apogee_estimates, apogee_estimate)
    
    plt.plot(timestamps, apogee_estimates, label="Estimated Apogee")
    plt.plot(timestamps, state_estimate, label="Current alitude")
    plt.show()
    
    