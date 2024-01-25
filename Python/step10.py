
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

########################################################
# Step 10: 2D Poisson Equation
########################################################

X = 50                       # Number of points along X-axis
Y = 50                       # Number of points along Y-axis
T = 100                      # Total number of time steps

xmin = 0   
xmax = 2
ymin = 0    
ymax = 1

dx = (xmax - xmin)/(X - 1.)
dy = (ymax - ymin)/(Y - 1.)

# Initialise spatial grid and initial condition
x = np.linspace(xmin, xmax, X)
y = np.linspace(ymin, ymax, Y)

p = np.zeros([Y, X])
b = np.zeros([Y, X])

b[int(Y / 4), int(X / 4)]  = 100
b[int(3 * Y / 4), int(3 * X / 4)] = -100
 
nX, nY = np.meshgrid(x, y)

# Setup the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surface = ax.plot_surface(nX, nY, p, cmap='viridis')
ax.set_title('2D Poisson Equation')
ax.set_zlim(-0.05, 0.05)

def animate(n):
    global p, surface
    pn = p.copy()
    for i in range(1, X-1):
        for j in range(1, Y-1):    
            p[i,j] = (dy * dy * (pn[i,j + 1] + pn[i,j - 1]) + dx * dx * (pn[i + 1,j] + pn[i - 1,j]) - b[i,j] * dy * dy * dx * dx) / (2 * (dy * dy + dx * dx))

    p[0, :] = 0
    p[:, 0] = 0

    p[Y-1, :] = 0
    p[:, X-1] = 0
   
    surface.remove()
    surface = ax.plot_surface(nX, nY, p, cmap='viridis')
    return surface,

ani = animation.FuncAnimation(fig, animate, frames=T, interval=70)

plt.show()