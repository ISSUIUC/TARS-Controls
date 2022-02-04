import numpy as np 
import pandas as pd
import src.atmosphere as atmosphere


def drag_from_csv(z, velocity_body, rasaero, Cd_list):
    # extracting the columns of interest
    mach_num = rasaero.mach.values; aoa = rasaero.alpha_deg.values; cd = rasaero.cd_power_off.values; protub = rasaero.protuberance.values
    # narrowing down the columns using mach number range (0.02 - 1.01)
    min_index = min(np.where(mach_num == 0.01)[0])
    max_index = min(np.where(mach_num == 1.01)[0])
    # re-make the lists using this range
    mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]
    mach = round(np.linalg.norm(velocity_body) / atmosphere.speed_sound(z),2)
    vx_b = velocity_body[0][0]
    vy_b = velocity_body[1][0]
    vz_b = velocity_body[2][0]
    alpha = abs(round(np.rad2deg(np.arctan2(vz_b,vx_b)), 0))

    mach_index_array = np.where(mach_num == mach)[0]; alpha_index_array = np.where(aoa == alpha)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))

    if len(intersection) == 0:
        return Cd_list[-1]

    index = intersection[0]
    return cd[index]

def drag_lookup_1dof(z,vel,rasaero,Cd_list):
    # extracting the columns of interest
    mach_num = rasaero.mach.values; aoa = rasaero.alpha_deg.values; cd = rasaero.cd_power_off.values; protub = rasaero.protuberance.values
    # narrowing down the columns using mach number range (0.01 - 1.01)
    min_index = min(np.where(mach_num == 0.01)[0])
    max_index = min(np.where(mach_num == 1.01)[0])
    # re-make the lists using this range
    mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]
    mach = round(vel / atmosphere.speed_sound(z),2)
 
    # Assuming vertical flight
    alpha = 0

    mach_index_array = np.where(mach_num == mach)[0]; alpha_index_array = np.where(aoa == alpha)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))

    if len(intersection) == 0:
        return Cd_list[-1]

    index = intersection[0]
    return cd[index]