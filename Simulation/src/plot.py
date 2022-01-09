import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# plots a trajectory based on an array of 3d points
def plot_3d(points):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = []
    x = []
    y = []
    for i in points:
        z.append(i[0])
        x.append(i[1])
        y.append(i[2])
    
    ax.scatter(x,y,z, c=[.75,0,0], alpha=1)
    

# plots a trajectory based on an array of 3d points 
# color maps the velocity based on a vector of velocity magnitudes
def plot_3d_vel(points, velocities):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = []
    x = []
    y = []
    for i in points:
        z.append(i[0])
        x.append(i[1])
        y.append(i[2])
    maxv = np.max(velocities)
    color = np.zeros((velocities.shape[0], 3))
    med = np.median(velocities)
    minv = np.min(velocities)

    # color points based on velocity
    for i in range(velocities.shape[0]):
        color[i, 0] = velocities[i]/maxv
        color[i, 1] = velocities[i]/maxv/4
        color[i, 2] = velocities[i]/maxv/4
    ax.scatter(x,y,z, c=color, alpha=1)

    # create color key
    max = mpatches.Patch(color=[1,0,0], label=maxv)
    m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
    m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
    m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
    min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
    ax.legend(handles=[max,m1,m2,m3,min])


# plots trajectory with velocity estimated based on a timestep
def plot_3d_est(points, timestep):
    velocities = np.zeros(points.shape[0])
    velocities[0] = abs(np.linalg.norm(points[1] - points[0]))/timestep
    for i in range(1, velocities.shape[0]):
        velocities[i] = abs(np.linalg.norm(points[i] - points[i-1]))/timestep
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    z = []
    x = []
    y = []
    for i in points:
        z.append(i[0])
        x.append(i[1])
        y.append(i[2])
    maxv = np.max(velocities)
    color = np.zeros((velocities.shape[0], 3))
    med = np.median(velocities)
    minv = np.min(velocities)

    # color points based on velocity
    for i in range(velocities.shape[0]):
        color[i, 0] = velocities[i]/maxv
        color[i, 1] = velocities[i]/maxv/4
        color[i, 2] = velocities[i]/maxv/4
    ax.scatter(x,y,z, c=color, alpha=1)

    # create color key
    max = mpatches.Patch(color=[1,0,0], label=maxv)
    m1 = mpatches.Patch(color=[1 - (.25*(minv/maxv)),0,0])
    m2 = mpatches.Patch(color=[1 - (.5*(minv/maxv)),0,0], label=med)
    m3 = mpatches.Patch(color=[1 - (.75*(minv/maxv)),0,0])
    min = mpatches.Patch(color=[minv/maxv,0,0], label= np.min(velocities))
    ax.legend(handles=[max,m1,m2,m3,min])