import matplotlib.pyplot as plt

kalman_dict = ...
color_arr = ["red","orange","yellow","green","blue", "indigo","violet","black","pink"]


plt.clf()
fig, axs = plt.subplots(2)
fig.suptitle('EKF Plots')

i = 0
for dim in ["x","y","z"]:
    for state in ["alt","vel","accel"]:
        axs[i%3].plot(kalman_dict[dim][state], color = color_arr[i]) 
        i += 1