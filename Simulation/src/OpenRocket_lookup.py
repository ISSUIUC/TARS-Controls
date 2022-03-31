import numpy as np 

def alt_vel_poly_fit(alt_des, ORK):
    time = ORK.Time; alt = ORK.Altitude; vel = ORK.Vertical_velocity

    min_index = min(np.where(time == 0)[0])
    max_index = min(np.where(time == 45.591)[0])
    time = time[min_index:max_index:1]; alt = alt[min_index:max_index:1]; vel = vel[min_index:max_index:1]
    poly = np.polyfit(alt, vel, 16)
    pfit = (np.poly1d(poly))(alt_des)
    return pfit