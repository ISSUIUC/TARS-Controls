import numpy as np
import src.rotation as rotation

# Rocket Properties
m0 = 21.364 #kg
r_CP = 219/100 #m
r_CG = 167.67/100 #m
I_rotational = 0.030245 #kg*m^2
I_longitudinal = 15.841 #kg*m^2
D = 0.1056132 #m

#Positions
x = 7857.1344 #m #initial alititude

PositionE = 116.21 #m
PositionN = 0.041193 #m
Lateral_Distance = 116.21 #m
Lateral_Direction = 0.02031 #degrees

#Velocities
vx = 309.25008 #m/s #initial veritcal velocity
lateral_velocity = 4.7643 #m/s

#Angular Rates
roll_rate = 1.73*10**(-6) #degrees/s #initial roll rate
pitch_rate = -0.01335 #degrees/s #initial pitch rate
yaw_rate = -0.0013 #degrees/s #initial yaw rate

#Accelerations
#gives vertical and lateral acceleration
ax = -16.442 #m/s^2 #initial vertical acceleration
lateral_accel = 0.19911 #m/s^2

#Orientation
yaw = 0.018162
pitch = 1.291
roll = 0

#Other Parameters
g = 9.8 #m/s^2
initial_alpha = 0.0043818 #deg
l_rocket = 3.02 #m
nose_ang = 0.069189 #nosecone angle (rad)

L_b = 2.2352 #Body Tube Length 
L_n = 0.762 #Nose cone length
T_f = 0.0029972 #Fin Thickness
L_m = 0.2032 #Length of fin from inner to root chord
n = 3 #Number of fins
A_fp = 0.011532235 #Fin Planform Artea
d_f = 0.08255 # Fin height

# Initialization
init_pos_f = np.array([[x],
                  [Lateral_Distance*np.cos(Lateral_Direction)],
                  [Lateral_Distance*np.sin(Lateral_Direction)]]) # X, Y, Z

init_or_f = np.array([[np.radians(yaw)],
                 [np.radians(pitch)],
                 [np.radians(roll)]]) # Yaw, Pitch, Roll

def init_vel_aligned_body(vx, lateral):
    # Get unit vector for xb axis
    uv_xb = np.array([[1],
                    [0],
                    [0]])

    # xb axis represented in the fixed frame
    R_fb = rotation.yaw(init_or_f[0][0]) @ rotation.pitch(init_or_f[1][0]) @ rotation.roll(init_or_f[2][0])
    uv_xb_f = R_fb @ uv_xb

    # Get divide lateral components up
    v_mag = (vx**2 + lateral**2)**(1/2)
    init_vel_f = uv_xb_f * v_mag
    return init_vel_f

init_vel_f = init_vel_aligned_body(vx, lateral_velocity)

init_angvel_f = np.array([[yaw_rate],
                        [pitch_rate],
                        [roll_rate]]) # Yaw rate, Pitch rate, Roll rate

