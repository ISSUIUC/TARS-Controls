from numpy import array


# values from openrocket mk3 file at mach 1
m0 = 21221  /1000 #kg
r_CP = 219  /100 #m
r_CG = 167.67  /100 #m
I_rotational = 0.030245 #kg*m^2
I_longitudinal = 15.841 #kg*m^2

#position
#gives altitude and east and north positions
x = 7712.5 #m #initial alititude

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
vx = 310.51 #m/s #initial veritcal velocity

lateral_velocity = 4.9938 #m/s
# vy = 
# vz = 
#v_0_f = np.array([vx,vy,vz]) # initial velocity in fixed frame

#angular velocity
roll_rate = 1.73*10**(-6) #degrees/s #initial roll rate
pitch_rate = -0.01335 #degrees/s #initial pitch rate
yaw_rate = -0.0013 #degrees/s #initial yaw rate
angular_vel_0_f = array([roll_rate,pitch_rate,yaw_rate]) #iniital angular velocities in fixed frame

#acceleration
#gives vertical and lateral acceleration
ax = -16.793 #m/s^2 #initital veritcal acceleration

lateral_accel = 0.19911 #m/s^2
#ay = 
#az = 
#a_0_f = np.array([ax,ay,az]) #initial acceleration in fixed frame

#angular acceleration
#got nothin