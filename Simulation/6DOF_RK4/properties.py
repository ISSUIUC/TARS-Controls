import math

# rocket mass w/out motor
mass = 1.0
# radius of rocket
r_r = 0.0254
# area w/out flaps
A = math.pi*r_r**2
# side profile area
A_s = 1.0
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6

# Temporary C_d constant 
#C_D's nuts in yo mouf
C_d = ...

# Motor Properties:
impulse = 17_907 # Ns
mass = 16.281 # Kg
delay = 0 # s
lookup_file = '../Lookup/Cesaroni_17907N2540-P_Trimmed.csv'