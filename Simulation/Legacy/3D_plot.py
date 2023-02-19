import src.plot_controls as hl
import matplotlib.pyplot as pyplt
import numpy as np 
# test script

plotted_hc = [[0,0,0], [1,1,0], [2,1.5,0], [3, 2, 0], [2, 2.5, 0]]
steps=400
step = .25
xvals = np.arange(-50, 50, step)
yvals = -(xvals**2) + 2500
z = np.zeros(steps)
plotted = np.zeros((steps, 3))
plotted[:, 1] = xvals
plotted[:, 0] = yvals
vel = abs(-2*xvals)
for i in range(vel.shape[0]):
    vel[i] += 20
    if vel[i] > 100:
        vel[i] = 100
print(vel)
# plot with velocity
hl.plot_3d_vel(plotted, vel)

hl.plot_3d_est(plotted, .25)

# plot without velocity
hl.plot_3d(plotted_hc)

pyplt.show()