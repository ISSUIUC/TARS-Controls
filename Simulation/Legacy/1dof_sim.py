from platform import mac_ver
import numpy as np
#import matplotlib.pyplot as plt
#from scipy import linalg
#from sympy import Matrix, false
#from mpl_toolkits import mplot3d
import pandas as pd
#from filterpy.common import Q_continuous_white_noise

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
#import src.plot_controls as plot
import src.rocket as rocket
#import src.rotation as rotation
import src.RASAero_lookup as rasaero
import src.altimeter as altimeter
import src.kalman_filter as kalman
import src.interpolation as interp
import src.propellant_mass as prop

#* ---------------------------- Frames we are using --------------------------- #
# Fixed Frame - fixed to the launch rail, not taking into account rotation of the Earth
# X points vertically
# Y is chosen through the right-hand rule
# Z points towards the launch rail

# Body-Fixed Frame - Origin at the center of mass of the rocket
# X points through the nose of the rocket
# Y is chosen through the right-hand rule
# Z points towards the launch lugs on the rocket

# Aerodynamic Frame (Body-Fixed) with Origin at the center-of-mass
# X aligned with aerodynamic velocity
# Y is chosen through the right-hand rule
# Z lies in plane that goes through the centerline of the rocket and the launch lug

#* ---------------------------- Naming Conventions ---------------------------- #
# pos = position(x,y,z)
# or = orientation(yaw,pitch,roll)
# vel = velocity(x,y,z)
# angvel = angular velocity(yaw,pitch,roll)
# accel = acceleration(x,y,z)
# angaccel = angular acceleration (yaw,pitch,roll)

#f = fixed frame
#b = body frame
#a = aerodynamic frame

# pos_f-> position in the fixed frame
# angvel_b-> angular velocity in the body frame
#* ------------------------------ Simulation Code ----------------------------- #


# Importing RasAero Package for Coeffiecient of Drag Lookup
#RASaero = pd.read_csv("Simulation/Lookup/RASAero.csv")
RASaero = pd.read_csv("C:/Users/kdcch/Documents/Github/TARS-Controls/Simulation/Lookup/RASAero_noAoA.csv")

# Import thrust curve csv
#thrust_csv = pd.read_csv("Simulation/Lookup/AeroTech_M2500T_Trimmed.csv")
# NOTE: if this csv is changed for a new motor, the burnout time in the if statements and mass constants will need to be changed
thrust_csv = pd.read_csv("C:/Users/kdcch/Documents/Github/TARS-Controls/Simulation/Lookup/AeroTech_M2500T_Trimmed.csv")

# Calculate moments of inertia and center of mass
#TODO: Move this into the simulation when simulating moving flaps -> Ixx changes
I, c_m, m = rocket.I_new(0,0)

s_dt = .003

# Initialize dictionary to store values at all time steps
ref_a_vels = []

dic = {
       "x":[],
       "x_noise":[],
       "vel": [],
       "accel": [],
       "CD": [],
       "Sref": [],
       "time_sim": [],
       "predict_alt": [],
       "predict_update_alt": []
       }

kalman_dic = {
        "alt": [],
        "vel": []
}

# Initial Values
# pos_f = constants.x
# pos_f_noisy = altimeter.alt_noise(constants.x)
# vel_f = constants.vx
# accel_f = constants.ax

pos_f = 0
pos_f_noisy = 0
vel_f = 0
accel_f = 0

#* Calculating the total mechanical energy, this is assumed to be a constant 
E_tot = m*9.81*pos_f + 0.5*m*vel_f**2

#* Kalman Filter Initialization
# Initialize states (x), measurement function (H), Covariance [P], White Noise [Q], Measurement Noise Function [R]
kalman.initialize(pos_f_noisy, vel_f, accel_f, s_dt)

# Time setup
start_time = 0.0
end_time = 45
total_steps = 10000
time = np.linspace(start_time,end_time,total_steps,endpoint=False)
dt = time[1] - time[0]

# Define max and min values for flap actuation
l_max = conversion.ft_to_m(1/12) # 1 inch actuation length
l_min = 0 # can't have negative actuation

# Define initial flap length at start of control time
l = 0
u = np.array([l])
for t in time:
    before_launch = True
    after_burnout = False
    print("altitude:",pos_f)
    print("velocity:",vel_f)


    # A-priori 
    kalman.priori(u)
    
    # Density varies with altitude
    rho = atmosphere.density(pos_f)

    # Calculate speed of sound at current altitude
    a = atmosphere.speed_sound(pos_f)

    # Calculate mach number at current time
    mach = vel_f/a

    # Total drag coefficient of airframe 
    #Cd_total = rasaero.drag_lookup_1dof(pos_f,vel_f,RASaero,dic["CD"])
    # Cd_total = 0
    Cd_total = interp.cd_interpolation(pos_f, vel_f, 0, l_max, l, RASaero, before_launch, after_burnout)

    #Approximation - use the area of a circle for reference area
    Sref_a = rocket.sref_approx(constants.D, l)
   
    if 0 <= t <= 0.019:
        thrust = 0
        #Calculate aerodynamic forces on the rocket and the acceleration in the aerodynamic frame
        F_a = 0
        accel_a = (thrust + F_a)/constants.m0 # Acceleration due to Aerodynamic Forces
        accel_f = 0 # Net Acceleration from Aerodynamic Forces + Gravity
    
    # add rocket forces for when motor is live during 0.082 < t < 4.264 s (for April Launch)
    # for final launch, 0.019 < t < 3.594
    elif 0.019 < t <= 3.594:
        before_launch = False
        #Calulate current mass of rocket at given time
        m_prop = prop.find_prop_mass_irec(t) # change function call to match which launch
        current_m = constants.m0 - m_prop

        thrust = interp.thrust_interp(t,thrust_csv)
        print("thrust:",thrust)

        #Calculate aerodynamic forces on the rocket and the acceleration in the aerodynamic frame
        F_a = -((rho* (vel_f**2) * Sref_a * Cd_total) / 2)
        accel_a = (thrust + F_a)/current_m # Acceleration due to Aerodynamic Forces
        print("accel:",accel_a)

        accel_f = accel_a - constants.g # Net Acceleration from Aerodynamic Forces + Gravity

    else:
        #Calculate aerodynamic forces on the rocket and the acceleration in the aerodynamic frame
        F_a = -((rho* (vel_f**2) * Sref_a * Cd_total) / 2)
        accel_a = F_a/constants.mf # Acceleration due to Aerodynamic Forces
        accel_f = accel_a - constants.g # Net Acceleration from Aerodynamic Forces + Gravity
        after_burnout = True

    #* End simulation if rocket reached apogee
    if (np.sign(vel_f) == -1):
        break

    #* Updating Values
    
    # Append Values to the Arrays
    dic["x"].append(float(pos_f))
    dic["x_noise"].append(float(pos_f_noisy))
    dic["vel"].append(float(vel_f))
    dic["accel"].append(float(accel_f))
    dic["CD"].append(float(Cd_total))
    dic["Sref"].append(float(Sref_a))
    dic["time_sim"].append(float(t))

    # Calculate new velocities and positions using current values
    pos_f = pos_f + (vel_f * dt) + (0.5 * (accel_f * (dt**2)))
    pos_f_noisy = altimeter.alt_noise(pos_f, mach)
    vel_f = vel_f + accel_f*dt
    
    # A-posteriori update

    kalman.update(pos_f_noisy, accel_f, Sref_a, rho)
    
    # kalman_dic["alt"].append(kalman.x_k[0][0])
    # kalman_dic["vel"].append(kalman.x_k[0][1])
    
#Print Apogee and total time taken
print("APOGEE (ft):", conversion.m_to_ft(max(dic["x"])))
print("Total Time Taken (s):", t)

#* --------------------------------- Plotting --------------------------------- #
#? Calculating the difference between the predicted altitude and actual altitude 
difference_pre = []
difference_post = []
for num in np.arange(0,len(dic["predict_alt"])):
    D = dic["predict_alt"][num] - max(dic["x"])
    d = dic["predict_update_alt"][num] - max(dic["x"])
    difference_pre.append(D)
    difference_post.append(d)



'''
# Position Measurements vs Kalman Filter Graph
plt.plot(dic["time_sim"], dic["x_noise"],label="Noisy Altitude Measurement",color="lightsteelblue",linestyle=":")
plt.plot(dic["time_sim"], dic["x"],label="True Altitude",color="royalblue")
plt.plot(dic["time_sim"], kalman.kalman_dic["alt"],label="Estimation",linestyle="--",color="tab:red")
plt.xlabel("Time (sec)")
plt.ylabel("Altitude (m)")
plt.legend()
plt.show()

# Velocity Measurements vs Kalman Filter Graph
plt.plot(dic["time_sim"], dic["vel"],label="True Velocity",color="royalblue")
plt.plot(dic["time_sim"], kalman.kalman_dic["vel"],label="Estimation",linestyle="--",color="tab:red")
plt.xlabel("Time (sec)")
plt.ylabel("Velocity (m/s)")
plt.legend()
plt.show()

# Acceleration Measurements vs Kalman Filter Graph
plt.plot(dic["time_sim"], dic["accel"],label="True Acceleration",color="royalblue")
plt.plot(dic["time_sim"], kalman.kalman_dic["accel"],label="Estimation",linestyle="--",color="tab:red")
plt.xlabel("Time (sec)")
plt.ylabel("Acceleration ($m/s^2$)")
plt.legend()
plt.show()

'''




# # Acceleration Plot
# plt.plot(time_flight,dic["accel"])
# plt.ylabel("Acceleration $(m / s^2)$")
# plt.xlabel("Time (s)")
# plt.show()

# # Velocity Plot
# plt.plot(time_flight,dic["vel"])
# plt.ylabel("Vertical Velocity $(m / s)$")
# plt.xlabel("Time (s)")
# plt.show()

# Altitude Plot
# plt.plot(time_flight,dic["x"])
# plt.plot(time_flight,dic["x_noise"])
# plt.plot(time_flight,kalman_dic["alt"])
# plt.ylabel("Altitude $(m)$")
# plt.xlabel("Time (s)")
# plt.show()

# # CD Plot
# plt.plot(time_flight,dic["CD"])
# plt.ylabel("CD")
# plt.xlabel("Time (s)")
# plt.show()