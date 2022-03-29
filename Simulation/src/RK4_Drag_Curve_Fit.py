import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
mach = np.arange(.01, .99, .01)
mach = np.round(mach, 2)
cdlist = np.empty(0)

RASaero = pd.read_csv("Simulation/Lookup/RASAero.csv")

mach_num = RASaero.mach.values; aoa = RASaero.alpha_deg.values; cd = RASaero.cd_power_off.values; protub = RASaero.protuberance.values
min_index = min(np.where(mach_num == 0.01)[0])
max_index = min(np.where(mach_num == 1.01)[0])
mach_num = mach_num[min_index:max_index:1]; aoa = aoa[min_index:max_index:1]; cd = cd[min_index:max_index:1]; protub = protub[min_index:max_index:1]
for i in mach:
    mach_index_array = np.where(mach_num == i)[0]; alpha_index_array = np.where(aoa == 0)[0]; protub_index_array = np.where(protub == 0)[0]
    mach_set = set(mach_index_array); alpha_set = set(alpha_index_array)
    intersection = list(mach_set.intersection(alpha_index_array))
    index = intersection[0]
    cdlist = np.append(cdlist, cd[index])

poly = np.polyfit(mach,cdlist, 18)
pfit = np.poly1d(poly)(mach)

diff = pfit - cdlist
maxDiff = np.amax(diff)

plt.plot(mach, cdlist, label = "lookup")
plt.plot(mach, pfit, label = "polyfit")
plt.xlabel("Mach Number")
plt.ylabel("Cd")
plt.legend()
plt.show()