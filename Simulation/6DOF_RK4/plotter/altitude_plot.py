import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.properties as prop
import properties.data_loader as dataloader


config = dataloader.config

import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["figure.autolayout"] = True
columns = ["time", "pos_x", "vel_x"]
df = pd.read_csv(os.path.join(os.path.dirname(__file__), config["meta"]["output_file"]), usecols=columns)
plt.plot(df.time, df.pos_x, label="altitude")
plt.plot(df.time, df.vel_x, label="vel_x")
plt.show()

