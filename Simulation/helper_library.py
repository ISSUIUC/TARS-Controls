import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix
from matplotlib import pyplot as plt
from statistics import mean
from mpl_toolkits import mplot3d
import matplotlib.pyplot as pyplt
import matplotlib.patches as mpatches
import numpy as np 
import math

from sympy.polys.polyroots import root_factors

class atmos:
    def density_func(z):
        #? Defining a few constants 
        #* temperature under standard condition (15 degrees C at sealevel) kelvin
        T_0 = 288.16 
        #* pressure under standard condition in (Pa)
        P_0 = 101325
        #* Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
        b = 0.0065
        #* gravitational constant 
        g = 9.81
        #* air density under standard condition (kg/m3)
        rho_0 = 1.225 
        #* ideal gas constant (J/kg K)
        R = 287.05

        #? returning the density at an altitude Z
        rho = rho_0*(1 - ((b*z)/T_0))**(g/(R*b))*(T_0/(T_0 - b*z))
        
        return rho
    def temp(z):
        #? Calculates temperature based on altitude 
        T_0 = 288.16
        #* Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
        b = 0.0065
        #* Current temperature based on altitute 
        T  = T_0 - b*z

        return T

    def viscosity(z):
        #? Takes in altitute (m) as input 
        u0 = 1.716e-5 
        Su = 111
        #* Reference temperature for Sutherland's Law
        T0 = 273.15
        #* temperature under standard condition (15 degrees C at sealevel) kelvin
        T_0 = 288.16
        #* Temperature lapse rate in k/m assuming temperature varies linearly based on altitude 
        b = 0.0065
        #* Current temperature based on altitute 
        T  = T_0 - b*z
        #* Calculates the dynamic viscosity
        miu = ((T/T0)**(1.5))*((T0 + Su)/(T + Su))*u0
        return miu
    
    def speed_sound(altitude):
        #* the first 11000 m of flight, the temperature varies linearly with altitude (k/m)
        dt_dh = 6.5e-3 
        #* specific heat ratio of air 
        r = 1.4 
        #* ideal gas constant of air R 
        R = 287.05
        #* temperature under standard condition (15 degrees C at sealevel) kelvin
        T_0 = 288.16 
        #* Current temperature based on altitute 
        T = T_0 - dt_dh*altitude
        #* Calculating speed of sound 
        a_H = np.sqrt(r*R*T)

        return a_H

class constants:
    #values from openrocket mk3 file at mach 1
    m0 = 21221  /1000 #kg
    r_CP = 219  /100 #m
    r_CG = 167.67  /100 #m
    I_rotational = 0.030245 #kg*m^2
    I_longitudinal = 15.841 #kg*m^2

    #position
    #gives altitude and east and north positions
    x = 7857.1344 #m #initial alititude

    PositionE = 116.21 #m
    PositionN = 0.041193 #m
    Lateral_Distance = 116.21 #m
    Lateral_Direction = 0.02031 #degrees
    #y = 
    #z = 
    #position_0_f = np.array([x,y,z]) #initial position in fixed frame

    #orientation
    #got nothing
    #phi = 
    #theta = 
    #phe = 
    #orientation_0_f = np.array([phi,theta,phe]) #initial orientation in fixed frame

    #velocity
    #gives vertical and laterall velocity
    vx = 309.25008 #m/s #initial veritcal velocity

    lateral_velocity = 4.7634144 #m/s
    # vy = 
    # vz = 
    #v_0_f = np.array([vx,vy,vz]) # initial velocity in fixed frame

    #angular velocity
    roll_rate = -2.54*10**(-8) # rad/s initial roll rate
    pitch_rate = 6.55*10**(-5) # rad/s initial pitch rate
    yaw_rate = -1.69*10**(-5) # rad/s initial yaw rate
    angular_vel_0_f = np.array([roll_rate,pitch_rate,yaw_rate]) #iniital angular velocities in fixed frame

    #acceleration
    #gives vertical and lateral acceleration
    ax = -16.793 #m/s^2 #initital veritcal acceleration

    lateral_accel = 0.19911 #m/s^2
    #ay = 
    #az = 
    #a_0_f = np.array([ax,ay,az]) #initial acceleration in fixed frame

    #angular acceleration
    #got nothin

class conversion:
    def ft_to_m(measurement): 
        return (measurement / 3.2808) 

    def c_to_k(measurement):
        return (measurement + 273.15)

    def deg_to_rad(measurement):
        return (measurement*np.pi/180)

    def rad_to_deg(measurement):
        return (measurement*180/np.pi)

class plot:
    
    # plots a trajectory based on an array of 3d points
    def plot_3d(points):
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = []
        x = []
        y = []
        for i in points:
            z.append(i[0])
            x.append(i[1])
            y.append(i[2])
        

        ax.scatter(x,y,z, c="tab:blue", alpha=1)
        ax.set_xlim3d(-60, 60); ax.set_ylim3d(-60, 60);

        return

    # plots a trajectory based on an array of 3d points 
    # color maps the velocity based on a vector of velocity magnitudes
    def plot_3d_vel(points, velocities):
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = []
        x = []
        y = []
        for i in points:
            z.append(i[0])
            x.append(i[1])
            y.append(i[2])
        maxv = np.max(velocities)
        color = np.zeros((velocities.shape[0], 3))
        med = np.median(velocities)
        minv = np.min(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
        m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
        m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
        min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
        ax.legend(handles=[max,m1,m2,m3,min])
        return

    # plots trajectory with velocity estimated based on a timestep
    def plot_3d_est(points, timestep):
        velocities = np.zeros(points.shape[0])
        velocities[0] = abs(np.linalg.norm(points[1] - points[0]))/timestep
        for i in range(1, velocities.shape[0]):
            velocities[i] = abs(np.linalg.norm(points[i] - points[i-1]))/timestep
        fig = pyplt.figure()
        ax = pyplt.axes(projection='3d')
        z = []
        x = []
        y = []
        for i in points:
            z.append(i[0])
            x.append(i[1])
            y.append(i[2])
        maxv = np.max(velocities)
        color = np.zeros((velocities.shape[0], 3))
        med = np.median(velocities)
        minv = np.min(velocities)

        # color points based on velocity
        for i in range(velocities.shape[0]):
            color[i, 0] = velocities[i]/maxv
            color[i, 1] = velocities[i]/maxv/4
            color[i, 2] = velocities[i]/maxv/4
        ax.scatter(x,y,z, c=color, alpha=1)

        # create color key
        max = mpatches.Patch(color=[1,0,0], label=maxv)
        m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
        m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
        m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
        min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
        ax.legend(handles=[max,m1,m2,m3,min])
        return

class inertia:
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

class rotation: 
    
    def yaw(phe):
        #* rotation about the z-axis, takes in angle phe in the function (rad)
        R_y = np.array([[np.cos(phe), -np.sin(phe), 0],[np.sin(phe), np.cos(phe), 0],[0, 0, 1]])
        return R_y

    def pitch(theta):
        #* rotation about the y-axis, takes in angle theta in the function (rad)
        R_p = np.array([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0],[-np.sin(theta), 0, np.cos(theta)]])
        return R_p
    
    def roll(phi):
        #* rotation about the x-axis, take sin angle phi in the function (rad)
        R_r = np.array([[1, 0, 0],[0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
        return R_r

    def body_aero(velocity_body):
        #? Takes in np.array velocity_body 
        #* rotation matrix that transforms the body frame to the aerodyanmic frame 
        #* takes in beta (side-slip angle), and alpha (angles of attack)
        #* needs the velocities in the body frame to calculate beta and alpha
        vx_b = velocity_body[0][0]
        vy_b = velocity_body[1][0]
        vz_b = velocity_body[2][0]

        alpha = np.arctan2(vz_b,vx_b)
        beta = np.arctan2(vy_b, np.sqrt(vx_b**2 + vz_b**2))
        R_ba = np.array([[np.cos(beta)*np.cos(alpha), np.sin(beta), np.cos(beta)*np.sin(alpha)], [-np.sin(beta)*np.cos(alpha), np.cos(beta), 
        -np.sin(beta)*np.sin(alpha)],[-np.sin(alpha), 0, np.cos(alpha)]])
        return R_ba
    
class sref:
    def sref_body(velocity_body,l,D):
        #? Takes in np.array velocity_body, the total length of rocket L, and the diameter of D as inputs 
        #* This function calculates the aerodynamic area of the rocket body
        #* Assuming the rocket body is a cylinder 
        #* This depends mainly on the sideslip angle beta 
        vx_b = velocity_body[0][0]
        vy_b = velocity_body[1][0]
        vz_b = velocity_body[2][0]

        beta = np.arctan2(vy_b, np.sqrt(vx_b**2 + vz_b**2))

        #* the length of the body sees by incoming air 
        l_effective = l*np.cos(beta)
        #* effective aerodynamic area 
        sref_b_eff = l_effective*D

        return l_effective, sref_b_eff

class coef_v1: 

    def friction_drag(z,l,D,velocity_body):
        #? Function takes in altitude "z", total length of the rocket "l", body tube diameter "D", and the body frame velocities 
        #? This function calculates the drag on the rocket due to viscous forces 
        #* This depends on the current Reynolds number and the critical Reynolds number 
        #* From the open rocket source http://openrocket.sourceforge.net/techdoc.pdf, and source no.4 listed in the
        #* reference, we find that the critical Reynolds number is arouond 5e5. 
        #? Import things to incorporate 
        #? 1) Viscosity Function 
        #? 2) Reynolds Number Equation 
        #? 3) Density-altitude function 

        rho = atmos.density_func(z)
        miu = atmos.viscosity(z)
        L = sref.sref_body(velocity_body,l,D)[0]
        v_mag = np.linalg.norm(velocity_body)
        #* Calculating the Reynold's number 
        Re = rho*v_mag*L/miu
        #* Critical Reynolds Number 
        # Re_c = 5e5
        Re_c = 1121800.2567142074
        #* Constant B in the equation in the paper: https://scholarcommons.scu.edu/cgi/viewcontent.cgi?article=1080&context=mech_senior
        B  = Re_c*((0.074/(Re**0.2)) - 1.328/(np.sqrt(Re)))
        #* Calculating Friction Coefficient 
        if Re <= Re_c:
            C_f = 1.328/(np.sqrt(Re))
        else: 
            C_f = (0.074/(Re**0.2)) - (B/Re)
        
        return C_f, Re

    def body_drag(l,L_b,L_n,d_b,C_f):
        #? l: total length of rocket
        #? L_b: length of body tube
        #? L_n: length of nose cone 
        #? d_b: diameter of body tube 
        #? d_d: diameter of base of rocket 
        
        Cd_body = (1 + (60/((l/d_b)**3)) + 0.0025*(L_b/d_b))*(2.7*(L_n/d_b) + 4*(L_b/d_b))*C_f
        return Cd_body
    
    def base_drag(d_b, d_d, Cd_body):
        #? d_b: diameter of body tube 
        #? d_d: diameter of base of rocket 

        Cd_base = 0.029*(d_b/d_d)**3/(np.sqrt(Cd_body))
        return Cd_base
    
    def fin_drag(T_f, L_m, n, A_fp, C_f, d_f):
        #? T_f: fin thickness
        #? L_m: true length of the fin from inner to outer edge 
        #? n: number of fins 
        #? A_fp: fin platform area 
        #? d_f: might be the width/height of the fins 
        #TODO: Double check this d_f value or find a new simpler equation 
        
        Cd_fin = 2*C_f*(1 + (T_f/L_m))*(4*n*A_fp)/(np.pi*d_f**2)

        return Cd_fin

    def interference_drag(T_f, L_m, n, A_fp, C_f, d_f, A_fe):
        #? T_f: fin thickness
        #? L_m: true length of the fin from inner to outer edge 
        #? n: number of fins 
        #? A_fp: fin platform area 
        #? d_f: might be the width/height of the fins 
        #? A_fe: exposed part of teh trapezoidal fin area 

        Cd_interfere = 2*C_f*(1 + (T_f/L_m))*(4*n*(A_fp - A_fe))/(np.pi*d_f**2)
        
        return Cd_interfere
        
class coef_v2: 

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

        rho = atmos.density_func(z)
        miu = atmos.viscosity(z)
        L = sref.sref_body(velocity_body,l,D)[0]
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
            C_f = 1.48*(10**(-2))
        elif Re > 10e4 and Re < Re_c:
            C_f = 1/((1.5*math.log(Re) - 5.6)**2)
        elif Re > Re_c:
            C_f = 0.032*((Rs/l)**0.2)
        
        #* Taking compressibility correction into account. Need correction factors 
        #* grabbing the speed of sound value
        a = atmos.speed_sound(z)
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
        a = atmos.speed_sound(z)
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
        a = atmos.speed_sound(z)
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

        Cd_total = coef_v2.friction_drag(z,l,D,velocity_body,A_ref) + coef_v2.pressure_base_drag(angle,z,velocity_body) + coef_v2.para_drag(z,velocity_body)

        return Cd_total

    def total_drag_scaled(z,l,D,velocity_body,A_ref,angle):

        A_nose = 0.29029
        A_fp = 0.011532235
        A_fins = 6*A_fp
        A_base = np.pi*((D/2)**2)
        A_para = 2*(1.1/100)*(1.6/100)

        Cd_nose, Cd_fins, Cd_base = coef_v2.pressure_base_drag(angle,z,velocity_body)
        
        Cd_scaled = (A_nose/A_ref)*Cd_nose + (A_fins/A_ref)*Cd_fins + (A_base/A_ref)*Cd_base + (A_para/A_ref)*coef_v2.para_drag(z,velocity_body) + coef_v2.friction_drag(z,l,D,velocity_body,A_ref)

        return Cd_scaled








