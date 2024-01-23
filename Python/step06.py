
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

########################################################
# Step 6: 2D Convection
########################################################

X = 31;                   # Number of points along X-axis
Y = 31;                   # Number of points along Y-axis
T = 40;                    # Total number of time steps

dx = 2. / (X - 1);         # Step size in the X direction
dy = 2. / (Y - 1);         # Step size in the Y direction
dt = 0.2 * dx;             # Time step size

# Initialise spatial grid and initial condition
x = np.linspace(0, 2, num=X)
y = np.linspace(0, 2, num=Y)

u = np.ones([Y, X])
v = np.ones([Y, X])

u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2 
v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2 

nX, nY = np.meshgrid(x, y)

# Setup the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surface = ax.plot_surface(nX, nY, u, cmap='viridis')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('2D Convection')

def animate(n):
    global u, v, surface
    un = u.copy()  
    vn = v.copy()  
    for i in range(1, X-1):
        for j in range(1, Y-1):    
            u[i,j] = un[i,j] - (un[i,j] * dt / dx * (un[i,j] - un[i-1,j])) - (vn[i,j] * (un[i,j] - un[i,j-1]) * dt / dy)
            v[i,j] = vn[i,j] - (un[i,j] * dt / dx * (vn[i,j] - vn[i-1,j])) - (vn[i,j] * (vn[i,j] - vn[i,j-1]) * dt / dy)
    u[0, :]  = 1
    u[:, 0]  = 1
    v[0, :]  = 1
    v[:, 0]  = 1
    
    u[-1, :] = 1
    u[:, -1] = 1
    v[-1, :] = 1
    v[:, -1] = 1
    

    surface.remove()
    surface = ax.plot_surface(nX, nY, u, cmap='viridis')
    return surface,

ani = animation.FuncAnimation(fig, animate, frames=T, interval=20)

plt.show()