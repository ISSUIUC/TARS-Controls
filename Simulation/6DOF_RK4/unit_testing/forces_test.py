'''
Creates a local copy of get_force and calls get_force to compare to known outputs in forces.py using rudimentary values. 
'''
import numpy as np
import sys
import os
import math
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
import dynamics.forces as forces 
import dynamics.rocket as rocket_model
import environment.atmosphere as atmosphere
import dynamics.motor as mot
import util.vectors as vct
import pandas as pd
import properties.properties as prop

# Defining neccessary variables  
cm = np.array([3.34-2.31, 0., 0.])
cp = np.array([3.34-2.71, 0., 0.])
rocket_dry_mass = 14.691
r_r = 0.0508
l = 3.34
A = math.pi*r_r**2
A_s = 2*r_r*l
max_ext_length = .0178
rasaero_file_location = os.path.join(os.path.dirname(__file__), prop.rasaero_lookup_file)
rasaero = pd.read_csv(rasaero_file_location)

# Create instances of neccessary objects to run get_force 
rocket = rocket_model.Rocket()
atm = atmosphere.Atmosphere() 
motor = rocket.motor
forces = forces.Forces(max_ext_length,cm,cp,A,A_s,rocket_dry_mass,motor,atm)
# Local "get_force copy" --> Designed to execute functions called within get_force
def get_force(x_state, flap_ext, time_stamp, density_noise=False) -> np.ndarray:
    alt = x_state.copy()[0,0]
    density = atm.get_density(alt, noise=density_noise, position=x_state[0])
    thrust = motor.get_thrust(time_stamp)
    wind_vector = atm.get_wind_vector(time_stamp)
    alpha = forces.get_alpha(x_state, wind_vector)
    test_get_alpha(alpha)
    drag = forces.aerodynamic_force(x_state, density, wind_vector, alpha, rasaero, thrust.dot(thrust) > 0, flap_ext) 
    values = forces.get_Ca_Cn_Cp(x_state, alpha, rasaero, thrust.dot(thrust) > 0, flap_ext)
    test_get_Ca_Cn_Cp(values)
    test_aerodynamic_force(drag)
    grav =  forces.gravitational_force(alt, time_stamp)
    test_gravitational_force(grav)
    force = vct.body_to_world(*x_state[2],thrust + drag) + grav
    moment = vct.body_to_world(*x_state[2], np.cross(-cm, thrust) + forces.aerodynamic_moment(drag))
    return np.array([force, moment]), alpha

# "get_force" from dynamics.forces
def test_get_force():
    #Initiallized Values (Random Values) (Working values of get_force as of 09/15/23)
    arr = np.array([[2,2,2],[2,2,2], [2,2,2]])
    answer = forces.get_force(arr, 0.015, 100.1)
    temp = get_force(arr, 0.015, 100.1)
    test_answer = np.array([[-152.8198926,-2.29827,-0.80334082],[0.74368931,-3.41707564,1.56364289]])

    #Testing output of get_force function
    # True --> the same 7-decimal points
    # False --> otherwise 
    test_array_component = True
    for i in range(0,np.size(answer[0],0)):
        for k in range(0, np.size(answer[0], 1)):
            boolean = (round(test_answer[i][k], 7) == round(answer[0][i][k], 7))
            if (boolean == False):
                test_array_component = False

    test_float_component = answer[1] == 0.9553166181245089
    response = (test_array_component and test_float_component)
    if response:
        return "PASSED --> [get_force]"
    else:
        return "FAILED --> [get_force]"

# "aerodynamic_force" from dynamics.forces
def test_aerodynamic_force(var):
    test_answer = np.array([0.06389041, -2.12639454, -8.94960151])
    passed = True 
    for i in range(0,3):
        if round(var[i], 7) != round(test_answer[i], 7):
            passed = False; 
    if passed:
        print("PASSED --> [aerodynamic_force]")
    else:
        print("FAILED --> [aerodynamic_force]")   

# "get_alpha" from dynamics.forces 
def test_get_alpha(var):
    answer = 0.9553166181245089
    if round(var, 7) == round(answer, 7):
        print("PASSED --> [get_alpha]")
    else:
        print("FAILED --> [get_alpha]")

# "gravitational_force" from dynamics.forces
def test_gravitational_force(var):
    answer = np.array([-143.94895119,0,0])
    response = True
    for i in range (0,3):
        if round(var[i], 7) != round(answer[i], 7):
            response = False
    if response:
        print("PASSED --> [gravitational_force]")
    else: 
        print("FAILED --> [gravitational_force]")

def test_get_Ca_Cn_Cp(var):
    answer0 = 2.59
    answer1 = 10.72
    answer2 = np.array([1.06356,0,0])
    response = True
    if (var[0] != answer0 or var[1] != answer1):
        response = False
        print("Fload Component Failed")
    if (not (np.allclose(answer2, var[2],0.00000001))):
        response = False
        print("Array Component Failed")
    if (response):
        print("PASSED --> [get_Ca_Cn_Cp]")
    else:
        print("FAILED --> [get_Ca_Cn_Cp]")




# Testing Input

print("--------------- \n Test Results | \n---------------")
print(test_get_force())

