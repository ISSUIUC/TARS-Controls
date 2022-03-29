import numpy as np
import src.rotation as rotation

# Mach 1 at Row 413

# Rocket Properties
m0 = 20.352 #kg
r_CP = 219/100 #m #TODO: Update
r_CG = 167.67/100 #m #TODO: Update
I_rotational = 0.030245 #kg*m^2 #TODO: Update
I_longitudinal = 15.841 #kg*m^2 #TODO: Update
D = 0.102 #m #! Check

#Positions
x = 8744.2 #m #initial alititude

PositionE = 116.21 #m #TODO: Update
PositionN = 0.041193 #m #TODO: Update
Lateral_Distance = 116.21 #m #TODO: Update
Lateral_Direction = 0.02031 #degrees #TODO: Update

#Velocities
vx = 305.58 #m/s #initial vertical velocity
lateral_velocity = 4.7643 #m/s # TODO: Update

#Angular Rates
roll_rate = 1.73*10**(-6) #degrees/s #initial roll rate #TODO: Update
pitch_rate = -0.01335 #degrees/s #initial pitch rate #TODO: Update
yaw_rate = -0.0013 #degrees/s #initial yaw rate #TODO: Update

#Accelerations
#gives vertical and lateral acceleration
ax = -15.605 #m/s^2 #initial vertical acceleration
lateral_accel = 0.19911 #m/s^2 #TODO: Update

#Orientation
yaw = 0.018162 #TODO: Update
pitch = 1.291 #TODO: Update
roll = 0 #TODO: Update

#Other Parameters
g = 9.8 #m/s^2
initial_alpha = 0.0043818 #deg # TODO: Update
l_rocket = 3.07 #m 
nose_ang = 0.069189 #nosecone angle (rad) #TODO: Update

#
#L_b = 2.2352 #Body Tube Length 
#L_n = 0.762 #Nose cone length
#T_f = 0.0029972 #Fin Thickness
#L_m = 0.2032 #Length of fin from inner to root chord
#n = 3 #Number of fins
#A_fp = 0.011532235 #Fin Planform Artea
#d_f = 0.08255 # Fin height

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

