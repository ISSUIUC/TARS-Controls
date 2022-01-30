from platform import mac_ver
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sympy import Matrix, false
from mpl_toolkits import mplot3d
import pandas as pd

#* Import Helper Function Library
import src.atmosphere as atmosphere
import src.constants as constants
import src.conversion as conversion
import src.plot_controls as plot
import src.rocket as rocket
import src.rotation as rotation

# NumPy arrays to store current states
pos_f = np.array([[constants.x],
                  [0],
                  [0]]) # X, Y, Z

or_f = np.array([[np.radians(-0.0030712)],
                 [np.radians(1.291)],
                 [0]]) # Yaw, Pitch, Roll

vel_f = np.array([[constants.vx],
                  [0],
                  [0]]) # Vx, Vy, Vz

# Get unit vector for xb axis
uv_xb = np.array([[1],
                  [0],
                  [0]])

# xb axis represented in the fixed frame
R_fb = rotation.yaw(or_f[0][0]) @ rotation.pitch(or_f[1][0]) @ rotation.roll(or_f[2][0])
uv_xb_f = R_fb @ uv_xb

# Get divide lateral components up


print(uv_xb_f)