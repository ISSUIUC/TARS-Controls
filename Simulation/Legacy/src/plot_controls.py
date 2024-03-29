import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# plots a trajectory based on a dictionary of 3d points
def plot_3d(dic, true_scale):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = (dic["x"])
    y = (dic["y"])
    x = (dic["z"])
    if (true_scale):
        ax.set_ylim([-np.max(z), np.max(z)])
        ax.set_xlim([-np.max(z), np.max(z)])

    ax.scatter(x,y,z, color = "tab:blue", alpha=1, label="Rocket Trajectory")
    ax.set_xlabel("X (m) ")
    plt.legend(fontsize = 15)
    # plt.ylim([-20,20]);plt.xlim([-20,20])
    
# plots a trajectory based on an array of 3d points 
# color maps the velocity based on a vector of velocity magnitudes
def plot_3d_vel(pos_dic, vel_dic, true_scale):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = (pos_dic["x"])
    y = (pos_dic["y"])
    x = (pos_dic["z"])
    velocities = []
    for i in range(len(vel_dic["Vx"])):
        velocities.append(np.sqrt(vel_dic["Vx"][i]**2 + vel_dic["Vy"][i]**2 + vel_dic["Vz"][i]**2))
    maxv = np.linalg.norm(np.max(velocities, axis=0))
    color = np.zeros((len(velocities), 3))
    med = np.linalg.norm(np.median(velocities, axis=0))
    minv = np.linalg.norm(np.min(velocities, axis=0))

    # color points based on velocity
    for i in range(len(velocities)):
        color[i, 0] = (np.linalg.norm(velocities[i])/maxv)
        if color[i, 0] > 1:
            color[i, 0] = 1
        color[i, 1] = (np.linalg.norm(velocities[i])/maxv)/4
        color[i, 2] = (np.linalg.norm(velocities[i])/maxv)/4
    if (true_scale):
        ax.set_ylim([-np.max(z), np.max(z)])
        ax.set_xlim([-np.max(z), np.max(z)])
    fig.suptitle("Velocity Based on User Input")
    ax.scatter(x,y,z, c=color, alpha=1)
    ax.set_xlabel("X (m) ")
    plt.legend(fontsize = 15)
    # create color key
    vrange = maxv-minv
    max = mpatches.Patch(color=[1,0,0], label=maxv)
    m1 = mpatches.Patch(color=[1 - (.25*(vrange))/maxv,0,0])
    m2 = mpatches.Patch(color=[1 - (.5*(vrange))/maxv,0,0], label=med)
    m3 = mpatches.Patch(color=[1 - (.75*(vrange))/maxv,0,0])
    min = mpatches.Patch(color=[minv/maxv,0,0], label=minv)
    ax.legend(handles=[max,m1,m2,m3,min])


# plots trajectory with velocity estimated based on a timestep
def plot_3d_est(pos_vals, timestep, true_scale):
    velocities = np.zeros(len(pos_vals["x"]))
    velocities[0] = abs(np.linalg.norm(pos_vals["y"] - pos_vals["x"]))/timestep
    for i in range(1, velocities.shape[0]):
        velocities[i] = abs(np.linalg.norm(pos_vals["x"][i] - pos_vals["x"][i-1]))/timestep
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = []
    x = []
    y = []
    for i in range(len(pos_vals["x"])):
        z.append(pos_vals["z"][i])
        x.append(pos_vals["x"][i])
        y.append(pos_vals["y"][i])
    maxv = np.max(velocities)
    color = np.zeros((velocities.shape[0], 3))
    med = np.median(velocities)
    minv = np.min(velocities)

    # color points based on velocity
    for i in range(velocities.shape[0]):
        color[i, 0] = velocities[i]/maxv
        color[i, 1] = velocities[i]/maxv/4
        color[i, 2] = velocities[i]/maxv/4
    if (true_scale):
        ax.set_ylim([-np.max(z), np.max(z)])
        ax.set_xlim([-np.max(z), np.max(z)])
    fig.suptitle("Estimated Velocity Based on Dt")
    ax.scatter(x,y,z, c=color, alpha=1)

    # create color key
    vrange = maxv-minv
    max = mpatches.Patch(color=[1,0,0], label=maxv)
    m1 = mpatches.Patch(color=[1 - (.25*(vrange))/maxv,0,0])
    m2 = mpatches.Patch(color=[1 - (.5*(vrange))/maxv,0,0], label=med)
    m3 = mpatches.Patch(color=[1 - (.75*(vrange))/maxv,0,0])
    min = mpatches.Patch(color=[minv/maxv,0,0], label=minv)
    ax.legend(handles=[max,m1,m2,m3,min])
    
def plot_accel_time(accel_vals, time_array):
    
    #Create figure and axes
    fig, axs = plt.subplots(2,sharex=True)
    fig.add_subplot(111,frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    
    #Storing values
    accel_x = []
    accel_y = []
    accel_z = []
    for val in range(len(accel_vals["Ax"])):
        accel_x.append(accel_vals["Ax"][val])
        accel_y.append(accel_vals["Ay"][val])
        accel_z.append(accel_vals["Az"][val])
        
    #Plotting
    axs[0].plot(time_array,accel_x,label="X Acceleration",color="maroon")
    axs[1].plot(time_array,accel_y,label="Y Acceleration")
    axs[1].plot(time_array,accel_z,label="Z Acceleration")
    axs[0].legend()
    axs[1].legend()
    fig.suptitle("Rocket Acceleration vs Time")
    plt.ylabel("Acceleration (m/s^2)")
    plt.xlabel("Time (s)")
    plt.show()
    
    
     