import numpy as np
import src.constants as constants
import src.atmosphere as atmosphere


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
    return(I,new_r_CG,new_m)

def sref(velocity_body,l,D):
    
    #TODO: Clean up comments  (atmosphere.py for reference)(atmosphere.py for reference)
    #? Takes in np.array velocity_body, the total length of rocket L, and the diameter of D as inputs 
    #* This function calculates the aerodynamic area of the rocket body
    #* Assuming the rocket body is a cylinder 
    #* This depends mainly on the sideslip angle beta 
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[1][0]
    vz_b = velocity_body[2][0]

    beta = np.arctan2(vy_b, np.sqrt(vx_b**2 + vz_b**2))
    alpha = np.arctan2(vz_b,vx_b)
    # true angle between body frame and aerodynamic frame x-axis
    gamma = np.sqrt(beta**2 + alpha**2)
    # effective aerodynamic reference area 
    sref_b_eff = l*np.sin(gamma)*D + np.pi*((D/2)**2)*np.cos(gamma)

    return sref_b_eff, beta

def sref_approx(D, flap_length):
    # Approximates reference area as a circle
    return np.pi*((D/2)**2) + 2*flap_length*0.0254

def friction_drag(z,l,D,velocity_body,A_ref):
    #? Function takes in altitude "z", total length of the rocket "l", body tube diameter "D", and the body frame velocities 
    #? This function calculates the drag on the rocket due to viscous forces 
    #* This depends on the current Reynolds number and the critical Reynolds number 
    #* From the open rocket source http://openrocket.sourceforge.net/techdoc.pdf, and source no.4 listed in the
    #* reference, we find that the critical Reynolds number is arouond 5e5. 
    #? Import things to incorporate 
    #? 1) Viscosity Function 
    #? 2) Reynolds Number Equation 
    #? 3) Density-altitude function 

    rho = atmosphere.density(z)
    miu = atmosphere.viscosity(z)
    L = sref_body(velocity_body,l,D)[0]
    v_mag = np.linalg.norm(velocity_body)
    #* Calculating the Reynold's number 
    Re = rho*v_mag*L/miu
    #* Critical Reynolds Number 
    # Re_c = 5e5
    #* Roughness Factor of an incorrectly sprayed aircraft surface Rs = 200 micro meter 
    Rs = 5e-6
    #* Critical Reynolds Number 
    Re_c = 51*(Rs/l)**(-1.039)

    if Re < 10e4:
        C_f = 1.48*10**(-2)
    elif Re > 10e4 and Re < Re_c:
        C_f = 1/((1.5*np.log(Re) - 5.6)**2)
    elif Re > Re_c:
        C_f = 0.032*((Rs/l)**0.2)
    
    #* Taking compressibility correction into account. Need correction factors 
    #* grabbing the speed of sound value
    a = atmosphere.speed_sound(z)
    #* Mach Number
    M = v_mag/a

    if M < 1:
        C_fc = C_f*(1 - 0.1*M**2)
    else:
        C_fc = C_f/((1 + 0.15*M**2)**(0.58))

    #* Calculating the Skin Friction Drag Coefficient 
    #? f_b: fineness ratio of rocket 
    f_b = l/D
    #? A_wet_body: wetted area, in contact with air, surface area of the nosecone + surface area of the body tube
    A_wet_cone = 0.29029
    A_wet_tube = np.pi*((D/2)**2)*2.2352
    A_wet_body = A_wet_tube + A_wet_cone
    #? A_wet_fins: wetted area of fins, twice the area of fins 
    #* A_fp: fin platform area
    A_fp = 0.011532235
    A_wet_fins = 6*A_fp
    #? A_ref: aerodynamic reference area of the rocket body

    #? t: thickness of the fins
    t = 0.0029972
    #? c: mean aerodynamic chord length 
    root_chord = 0.2032
    tip_chord = 0.0762
    taper_ratio = tip_chord/root_chord
    c = root_chord*(2/3)*(1 + taper_ratio + taper_ratio**2)/(1 + taper_ratio)
    Cd_f = C_fc * ((1 + (1/(2*f_b)))*A_wet_body + (1 + 2*t/c)*A_wet_fins)/A_ref

    return Cd_f

def pressure_base_drag(angle,z,velocity_body):
    #? First calculate the pressure drag due to the nosecone, nosecone joint angle is 0.069189
    Cd_pres_nose = 0.8*(np.sin(angle)**2)
    #? Next calculate the fin pressure drag, assuming square profile, fin pressure will equal base drag 
    #* Calculate Mach number 
    v_mag = np.linalg.norm(velocity_body)
    #* Speed of sound 
    a = atmosphere.speed_sound(z)
    #* Mach number 
    M = v_mag/a

    if M < 1:
        Cd_base = 0.12 + 0.13*M**2
    else:
        Cd_base = 0.25/M

    Cd_fin = Cd_base

    Cd_pres_base = Cd_pres_nose + Cd_fin + Cd_base

    return Cd_pres_nose, Cd_fin, Cd_base

def para_drag(z,velocity_body):
    #* Calculate Mach number 
    v_mag = np.linalg.norm(velocity_body)
    #* Speed of sound 
    a = atmosphere.speed_sound(z)
    #* Mach number 
    M = v_mag/a
    #? Stagpressure coefficient Calculation 
    if M < 1:
        dynamic_pres_ratio = 1 + M**2/4 + M**4/40
    else:
        dynamic_pres_ratio = 1.84 - (0.76/(M**2)) + (0.166/(M**4)) + (0.035/(M**6))

    Cd_stag = 0.85 * dynamic_pres_ratio
    #? Parasitic drag
    l_over_d = 1.1/1.6
    Cd_para = max(1.3 - 0.3*l_over_d,1)*Cd_stag

    return Cd_para

def total_drag(z,l,D,velocity_body,A_ref,angle):
    
    Cd_total = friction_drag(z,l,D,velocity_body,A_ref) + pressure_base_drag(angle,z,velocity_body) + para_drag(z,velocity_body)
    
    return Cd_total

def total_drag_scaled(z,l,D,velocity_body,A_ref,angle):

    A_nose = 0.29029
    A_fp = 0.011532235
    A_fins = 6*A_fp
    A_base = np.pi*((D/2)**2)
    A_para = 2*(1.1/100)*(1.6/100)

    Cd_nose, Cd_fins, Cd_base = pressure_base_drag(angle,z,velocity_body)
    
    Cd_scaled = (A_nose/A_ref)*Cd_nose + (A_fins/A_ref)*Cd_fins + (A_base/A_ref)*Cd_base + (A_para/A_ref)*para_drag(z,velocity_body) + friction_drag(z,l,D,velocity_body,A_ref)

    return Cd_scaled

