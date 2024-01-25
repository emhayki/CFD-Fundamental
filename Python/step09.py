
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

########################################################
# Step 9: Laplace equation
########################################################

X = 31;                       # Number of points along X-axis
Y = 31;                       # Number of points along Y-axis
T = 500;                      # Total number of time steps

dx = 2. / (X - 1);            # Step size in the X direction
dy = 2. / (Y - 1);            # Step size in the Y direction

# Initialise spatial grid and initial condition
x = np.linspace(0, 2, num=X)
y = np.linspace(0, 2, num=Y)

p = np.zeros([Y, X])
p[-1, :] = y

nX, nY = np.meshgrid(x, y)

# Setup the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surface = ax.plot_surface(nX, nY, p, cmap='viridis')
ax.set_title('Laplace equation')
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.view_init(30, 225)

def animate(n):
    global p, surface
    pn = p.copy()
    for i in range(1, X-1):
        for j in range(1, Y-1):    
            p[i,j] = (dy * dy * (pn[i+1,j]  + pn[i-1,j] ) + dx * dx * (pn[i,j+1]  + pn[i,j-1] ))/(2 * (dy * dy  + dx * dx))

    p[:,  0] = 0  
    p[:, -1] = y
    p[0,  :] = p[1, :] 
    p[-1, :] = p[-2, :]
   
    surface.remove()
    surface = ax.plot_surface(nX, nY, p, cmap='viridis')
    return surface,

ani = animation.FuncAnimation(fig, animate, frames=T, interval=70)

plt.show()