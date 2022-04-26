import numpy as np
import src.rotation as rotation

# Rocket Properties
m0_IREC = 30.544 # kg, for IREC, full mass before burnout
m0_April = 23.782 # kg, for April, full mass before motor burn
mf_IREC = 21.364 #kg, for IREC launch, after burnout
mf_April = 19.0586 # kg, for April, after burnout

r_CP = 219/100 #m
r_CG = 167.67/100 #m
I_rotational = 0.030245 #kg*m^2
I_longitudinal = 15.841 #kg*m^2
D = 0.1056132 #m


April_apogee_goal = 14000  # ft
IREC_apogee_goal = 30000   # ft

April_m = 19.0586
IREC_m = 21.1066

# Accelerations at start of burnout (not needed)
# April_accel = -31.593 #m/s^2
# IREC_accel = -30 #! FIX THIS

# Times at which thrust starts and ends
thrust_start_April = 0.082
thrust_start_IREC = 0.019
thrust_end_April = 4.264
thrust_end_IREC = 3.594

#Positions
x = 0 #m #initial alititude (burnout)

PositionE = 116.21 #m #TODO: Update
PositionN = 0.041193 #m #TODO: Update
Lateral_Distance = 116.21 #m #TODO: Update
Lateral_Direction = 0.02031 #degrees #TODO: Update

#Velocities
vx = 0
lateral_velocity = 4.7643 #m/s # TODO: Update

#Angular Rates
roll_rate = 1.73*10**(-6) #degrees/s #initial roll rate #TODO: Update
pitch_rate = -0.01335 #degrees/s #initial pitch rate #TODO: Update
yaw_rate = -0.0013 #degrees/s #initial yaw rate #TODO: Update

#Accelerations
#gives vertical and lateral acceleration
ax = 0 #m/s^2 #initial vertical acceleration
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
flap_width = 0.0351 #m
max_flap_length = .816 #in
gamma = 1.4 # specific heat ratio

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

