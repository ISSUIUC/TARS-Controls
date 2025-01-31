import math
import numpy as np

# Config file
sim_config = "../properties/configs/stargazer14_flightlike.yaml"
tiltlock_startkey = "    - sustainer:"
tiltlock_endkey = "desired_apogee: 1824"


### Calculated constants
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6


