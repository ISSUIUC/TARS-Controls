import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import pandas as pd
import math

import dynamics.forces as forces
import util.vectors as vct
import properties.properties as prop
import dynamics.rocket as rocket_model
class Apogee: 
    '''Apogee Estimator Class
    
    Args:
        state (np.array): 2D array of state vectors --> ([pos], [vel], [acc])
        dt (float): time step for RK4
        a (float): lower bound of domain (inclusive)
        b (float): upper bound of domain (inclusive)
        n (int): number of points to use in approximation
        atm (Atmosphere): atmosphere object
        simulator_config: Simulator configuration file
    '''

    atm = None
    rasaero_file_location = ""
    rasaero = None
    rocket = None
    stage_config = None
    
    def __init__(self, state, dt, a, b, n, atm, stage_config):
        self.rasaero_file_location = os.path.join(os.path.dirname(__file__), stage_config['rocket_body']['rasaero_lookup_file'])
        self.rasaero = pd.read_csv(self.rasaero_file_location)
        self.state = state[:3][0].copy()
        self.dt = dt
        self.flap_ext = 0.
        self.c, self.x_interpolate = self.calc_spline_coefficients(a, b, n)
        self.n = n
        self.atm = atm
        self.stage_config = stage_config
        self.rocket = rocket_model.Rocket(self.stage_config)
    
    def set_params(self, state):
        '''Reset the state vector to the current state of the rocket
        
        Args:
            state (np.array): 2D array of state vectors --> ([pos], [vel], [acc])
        '''
        self.state = state.copy()

    def RK4(self):
        """Runge-Kutta 4th order method for state propogation
        """
        k1_v = self.state[2]
        k2_v = self.state[2] + k1_v*(self.dt/2)
        k3_v = self.state[2] + k2_v*(self.dt/2)
        k4_v = self.state[2] + k3_v*(self.dt)


        v = (self.state[1] + (1/6)*(k1_v+(2*k2_v)+(2*k3_v)+k4_v)*self.dt)
    
        k1_p = self.state[1]
        k2_p = self.state[1] + k1_p*(self.dt/2)
        k3_p = self.state[1] + k2_p*(self.dt/2)
        k4_p = self.state[1] + k3_p*(self.dt)

        p = (self.state[0] + (1/6)*(k1_p+(2*k2_p)+(2*k3_p)+k4_p)*self.dt)

        a = self.get_accel()

        self.state = np.array([p,v,a])
    
    def step_p(self, y0, y1, dt):
        """Calculates slope of position over given delta t for state propogation
        
        Args:
            y0 (float): initial position
            y1 (float): final position
            dt (float): time step
            
        Returns:
            (float): slope of position over given delta t
        """
        return (y1-y0)/dt
    
    def get_Ca(self) -> float:
        # TODO: account for area change of flaps in C_a calculation
        '''References lookup table to find C_a based on flap extension

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

        Ca = self.approximate_cubic_spline(self.c, self.x_interpolate, mach)
        return Ca


    def get_accel(self) -> np.ndarray:
        '''Calculates net force felt by rocket while accounting for thrust, drag, gravity, wind
        
        Returns:
            (np.array): 2D array of forces and moments --> ([Fx, Fy, Fz], [Mx, My, Mz])
        '''
        # TODO: Add random disturbances
        # print("State: ", x_state)
        C_a = self.get_Ca()
        alt = self.state[0]
        density = self.atm.get_density(alt)
        drag = -0.5*(self.state[1]**2 * C_a*density*self.rocket.A)
        grav = self.gravitational_force(alt)
        force =  drag + grav
        return force/self.rocket.rocket_dry_mass

    def gravitational_force(self, altitude) -> np.ndarray:
        '''Calculates gravitational force acting on rocket based on altitude
        Relevant Equations:
            F = GMm/r^2

        Args:
            altitude (float): current altitude of rocket
        
        Returns:
            (np.array): vector of gravitational forces on each axis [1x3]
        '''
        return -(prop.G*prop.m_e*self.rocket.rocket_dry_mass)/((prop.r_e+altitude)**2)

    def step_v(self):
        
        '''Calculates slope of v over given delta t for state propogation
        
        Returns:
            (np.array): rate of change of velocity (acceleration) in form of state vector
        '''
        return self.get_accel()    
    
    def predict_apogee(self, current_state):
        '''Runs RK4 to predict apogee of rocket given current state
        
        Args:
            current_state (np.array): current state of rocket in form of state vector
            
        Returns:
            self.state[0] (float): predicted apogee of rocket
        '''
        self.set_params(current_state.copy())
        while (self.state[1] > 0):
            self.RK4()
        return self.state[0]
    
    ########################################
    # C_a function approximation functions #
    ########################################
    def f_true(self, x):
        x = round(x,2)
        return self.rasaero[self.rasaero['Mach Number'] == x]['CA Power-Off'].values[0]
    
    def calc_spline_coefficients(self, a, b, n):
        """Returns function approximations and error for natural cubic spline interpolation.

        This function assumes an f_true(x) function is globally available for calculating the true function value at x.
        
        Args:
            a (float): lower bound of interpolation interval
            b (float): upper bound of interpolation interval
            n (int): number of interpolation points
            
        Returns:
            (np.array): 4xn array of spline coefficients
        """
        
        # get interpolation points (uniform) and delta x
        x_interpolate = np.linspace(a, b, n+1).tolist()
        dx = x_interpolate[1] - x_interpolate[0]

        # get A matrix and f vector
        A = np.zeros((4 * n, 4 * n))
        f = np.zeros(4 * n)
        for i in range(n): # loop through each local spline (i)

            # get first row/column index associated with the ith spline
            ind = i * 4
            
            # update values from condition (1), (5)
            A[ind, ind] = dx**2/6
            A[ind, ind + 1] = 0
            A[ind, ind + 2] = x_interpolate[i]
            A[ind, ind + 3] = 1

            f[ind] = self.f_true(x_interpolate[i])

            # update values from condition (2), (6)
            A[ind + 1, ind] =  0
            A[ind + 1, ind + 1] = dx**2/6
            A[ind + 1, ind + 2] = x_interpolate[i+1]
            A[ind + 1, ind + 3] = 1

            f[ind + 1] = self.f_true(x_interpolate[i+1])
            
            if i == n - 1:
                # update values from "extra" condition (7)
                A[ind + 2, 0] = 1

                # update values from "extra" condition (8)
                A[ind + 3, ind + 1] = 1
            else:
                # update values from condition (3)
                A[ind + 2, ind] = 0
                A[ind + 2, ind + 1] = dx/2
                A[ind + 2, ind + 2] = 1
                A[ind + 2, ind + 3] = 0
                A[ind + 2, ind + 4] = dx/2
                A[ind + 2, ind + 5] = 0
                A[ind + 2, ind + 6] = -1
                A[ind + 2, ind + 7] = 0

                # update values from condition (4)
                A[ind + 3, ind] = 0
                A[ind + 3, ind + 1] = 1
                A[ind + 3, ind + 4] = -1

            # update values from conditions (3), (4) and "extra" conditions (7)-(8)
            f[ind + 2] = 0 # f_true(0)
            f[ind + 3] = 0 # f_true(0)
            

        # solve matrix system
        c = np.linalg.solve(A, f)
        return c, x_interpolate
    
    def approximate_cubic_spline(self, c, x_interpolate, x):
        """
        Returns the approximated function value for x using the cubic spline coefficients c and the interpolation points x_interpolate.
        
        Args:
            c (array_like): Vector of cubic spline coefficients
            x_interpolate (array_like): List of interpolation points used
            x (float): Value to evaluate the approximated function at
        
        Returns:
            fa_val (float_like): Approximated function value at x
        """
        # get spline index i for x
        if x == x_interpolate[-1]:
            i = len(x_interpolate) - 1   # index for last interpolation point
        elif x in x_interpolate:
            i = x_interpolate.index(x)  # index for interpolation points
        else:
            i = math.floor(x*(self.n/3)) # index for test points between interpolation points

        # get first row index associated with the ith spline in c
        ind = i * 4
        
        # get spline i output ("\"" is a line continuation character)
        fa_val = c[ind]/(6*(x_interpolate[i]-x_interpolate[i+1]))*(x-x_interpolate[i+1])**3 + \
                    c[ind+1]/(6*(x_interpolate[i+1]-x_interpolate[i]))*(x-x_interpolate[i])**3 + \
                    c[ind+2]*x + c[ind+3]
        return fa_val

if __name__=='__main__':
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    data = pd.read_csv(os.path.join(os.path.dirname(__file__), "LookUp", "flight_computer_20221029.csv"))
    state_estimate_x = data["state_est_x"]
    state_estimate_vx = data["state_est_vx"]
    state_estimate_ax = data["state_est_ax"]
    state_estimate = np.array([state_estimate_x, state_estimate_vx, state_estimate_ax])
    timestep = 0.1
    ax = data["ax"]
    rocket_estimated_apogee = data["state_est_apo"]
    timestamps = data["timestamp_ms"]
    max_timestamp = 5000
    timestamps = timestamps[:max_timestamp:int(timestep*100)]

    apogee_estimator = Apogee(state_estimate, timestep, 0.01, 3, 300)
    apogee_estimates = np.array([])
    for estimate in tqdm(state_estimate.T[:max_timestamp:int(timestep*100)]):
        apogee_estimate = apogee_estimator.predict_apogee(estimate)
        apogee_estimates = np.append(apogee_estimates, apogee_estimate)
    
    plt.plot(timestamps, apogee_estimates, label="Estimated Apogee")
    plt.plot(timestamps, state_estimate_x[:max_timestamp:int(timestep*100)], label="Current alitude")
    plt.plot(timestamps, rocket_estimated_apogee[:max_timestamp:int(timestep*100)], label="Rocket Estimated Apogee")
    plt.legend()
    
    plt.figure()
    true_apo = np.amax(state_estimate_x)
    apogee_estimate_error_absolute = (true_apo * np.ones_like(state_estimate_x[:max_timestamp:int(1/timestep)]) - apogee_estimates)
    apogee_estimate_error_percent = apogee_estimate_error_absolute/true_apo
    plt.plot(timestamps, apogee_estimate_error_absolute, label="Absolute Error")
    plt.legend()
    plt.figure()
    plt.plot(timestamps, apogee_estimate_error_percent, label="Percent Error")
    plt.legend()
    plt.show()
    
    
