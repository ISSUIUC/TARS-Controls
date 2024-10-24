import numpy as np 
import pandas as pd
import src.atmosphere as atmosphere
import warnings
warnings.simplefilter('ignore', np.RankWarning)


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

def drag_lookup_1dof(z,vel,rasaero,Cd_list, input, before_launch, before_burnout):
    
    # No drag on ground
    if (before_launch):
        return 0
    
    # Choose CD list based on power
    cd = rasaero.CD_Power_Off.values
    if (before_burnout):
        cd = rasaero.CD_Power_On.values
            
    # extracting the columns of interest
    mach_num = rasaero.mach.values; aoa = rasaero.alpha.values; cd = rasaero.CD_Power_Off.values; protub = rasaero.protuberance.values
    # narrowing down the columns using mach number range (0.01 - 2.50)
    min_index = min(np.where(mach_num == 0.01)[0])
    max_index = min(np.where(mach_num == 2.50)[0])
    # re-make the lists using this range
    mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]
    mach = round(vel / atmosphere.speed_sound(z),2)
 
    # Assuming vertical flight
    alpha = 0

    mach_index_array = np.where(mach_num == mach)[0]; alpha_index_array = np.where(aoa == alpha)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))
    
    # Correction     

    # No number found, use the previous value
    if len(intersection) == 0:
        if (mach < 0.01):
            return 0
        return Cd_list[-1]

    # index = intersection[0]
    index1 = intersection[0]; index2 = intersection[1]
    
    Cd = cd[index1] + input*39.3701*(cd[index2] - cd[index1])
    return Cd

def drag_lookup_curve_fit_poly(with_thrust, RASaero):
    mach = np.arange(.01, 2.49, .01)
    mach = np.round(mach, 2)
    cdlist = np.empty(0)

    RASaero = pd.read_csv("Simulation/Lookup/RASAero_Intrepid_5800_mk6.csv")

    mach_num = RASaero.mach.values; aoa = RASaero.alpha.values; cd = RASaero.CD_Power_Off.values; protub = RASaero.protuberance.values
    if (with_thrust):
        cd = RASaero.CD_Power_On.values
        
    min_index = min(np.where(mach_num == 0.01)[0])
    max_index = min(np.where(mach_num == 2.50)[0])
    mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]
    for i in mach:
        mach_index_array = np.where(mach_num == i)[0]; alpha_index_array = np.where(aoa == 0)[0]; protub_index_array = np.where(protub == 0)[0]
        mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
        intersection = list(mach_set.intersection(alpha_index_array))
        index = intersection[0]
        cdlist = np.append(cdlist, cd[index])

    poly = np.polyfit(mach,cdlist, 150)
    return poly