from numpy import array
# values from openrocket mk3 file at mach 1
m0 = 21.364 #kg
r_CP = 219  /100 #m
r_CG = 167.67  /100 #m
I_rotational = 0.030245 #kg*m^2
I_longitudinal = 15.841 #kg*m^2

#position
x = 7857.1344 #m #initial alititude

PositionE = 116.21 #m
PositionN = 0.041193 #m
Lateral_Distance = 116.21 #m
Lateral_Direction = 0.02031 #degrees

vx = 309.25008 #m/s #initial veritcal velocity

lateral_velocity = 4.763 #m/s


#angular velocity
roll_rate = 1.73*10**(-6) #degrees/s #initial roll rate
pitch_rate = -0.01335 #degrees/s #initial pitch rate
yaw_rate = -0.0013 #degrees/s #initial yaw rate
angular_vel_0_f = array([roll_rate,pitch_rate,yaw_rate]) #iniital angular velocities in fixed frame

#acceleration
#gives vertical and lateral acceleration
ax = -16.442 #m/s^2 #initital veritcal acceleration

lateral_accel = 0.19911 #m/s^2
