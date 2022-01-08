import numpy as np
import constants
import atmosphere


def I_new(m,r_m):
    
    # m is added mass (in grams) to nosecone 
    # r_m is the distance (in cm) from the tip of the nosecone to the CG of added mass (can be estimated) 

    # m = input("dry mass addition(kg):");m = float(m)
    # r_m = input("dry mass location(m):");r_m = float(r_m)

    #constants obtained from open rocket data
    r_CG_0 = constants.r_CG 
    m0 = constants.m0
    Ixx_0 = constants.I_rotational 
    Iyy_0 = constants.I_longitudinal
    Izz_0 = Iyy_0
    
    #finding new center of mass location
    new_m = m0 + m
    new_r_CG = (m0*r_CG_0 + m*r_m) / new_m
    d1 = new_r_CG - r_m #distance between new CG and the location where the mass is added
    d2 = r_CG_0 - new_r_CG #distance between new and old CGs
    
    #calculating new moments of inertia around new CG 
    Ixx_new = Ixx_0
    Iyy_new = Iyy_0 + m*d1**2 + m0*d2**2
    Izz_new = Iyy_new
    I = np.diag(np.array([Ixx_new,Iyy_new,Izz_new]))
    return(I,new_r_CG)

def sref(velocity_body,l,D):
    
    #? Takes in np.array velocity_body, the total length of rocket L, and the diameter of D as inputs 
    #* This function calculates the aerodynamic area of the rocket body
    #* Assuming the rocket body is a cylinder 
    #* This depends mainly on the sideslip angle beta 
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[0][1]
    vz_b = velocity_body[0][2]

    beta = np.arctan2(vy_b, np.sqrt(vx_b**2 + vz_b**2))

    #* the length of the body sees by incoming air 
    l_effective = l*np.cos(beta)
    #* effective aerodynamic area 
    sref_b_eff = l_effective*D

    return l_effective, sref_b_eff

def friction(z,l,D,velocity_body):
    
    # TODO: Clean up comments
    #? Function takes in altitude "z", 
    #? This function calculates the drag on the rocket due to viscous forces 
    #* This depends on the current Reynolds number and the critical Reynolds number 
    #* From the open rocket source http://openrocket.sourceforge.net/techdoc.pdf, and source no.4 listed in the
    #* reference, we find that the critical Reynolds number is arouond 5e5. 
    #? Import things to incorporate 
    #? 1) Viscosity Function 
    #? 2) Reynolds Number Equation 
    #? 3) Density-altitude function 

    rho = atmosphere.density(z)
    mu = atmosphere.viscosity(z)
    L = sref.sref_body(velocity_body,l,D)[0]
    v_mag = np.linalg.norm(velocity_body)
    
    #* Calculating the Reynold's number 
    Re = rho*v_mag*L/mu
    #* Critical Reynolds Number 
    Re_c = 5e5
    #* Constant B in the equation in the paper: https://scholarcommons.scu.edu/cgi/viewcontent.cgi?article=1080&context=mech_senior
    B  = Re_c*((0.074/(Re**0.2)) - 1.328/(np.sqrt(Re)))
    #* Calculating Friction Coefficient 
    if Re <= Re_c:
        C_f = 1.328/(np.sqrt(Re))
    else: 
        C_f = 0.074/(Re**0.2) - B/Re
    
    return C_f