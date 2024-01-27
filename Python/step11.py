
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

########################################################
# Step 11: Cavity Flow 
########################################################

X = 41
Y = 41
T = 500

dx = 2 / (X - 1)
dy = 2 / (Y - 1)

rho = 1
nu  = 0.1
dt  = 0.001

p, b = np.zeros((Y, X)), np.zeros((Y, X))
u, v = np.zeros((Y, X)), np.zeros((Y, X))

x = np.linspace(0, 2, X)
y = np.linspace(0, 2, Y)
nX, nY = np.meshgrid(x, y)

# Lid-driven cavity condition
u[-1, :] = 1 

# Setup the plot
fig, ax = plt.subplots()
speed   = np.sqrt(u**2 + v**2)
contour = ax.contourf(nX, nY, p, cmap='jet')

streamplot = ax.streamplot(nX, nY, u, v, color=speed, linewidth=0.5, cmap='jet')
cbar       = fig.colorbar(streamplot.lines)
ax.set_ylim([0, 2.25])

# Pressure Poisson equation
def pressure_poisson(p, dx, dy, b):
    for _ in range(50):
        pn = p.copy()
        for i in range(1, Y - 1):
            for j in range(1, X - 1):
                p[i, j] = (((pn[i, j+1] + pn[i, j-1]) * dy**2 +
                            (pn[i+1, j] + pn[i-1, j]) * dx**2) /
                           (2 * (dx**2 + dy**2)) -
                           dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[i, j])

        # Boundary conditions
        p[:, -1], p[0, :], p[:, 0] = p[:, -2], p[1, :], p[:, 1]
        p[-1, :] = 0
    return p

def animate(n):
    global u, v, p

    for i in range(1, Y - 1):
        for j in range(1, X - 1):
            b[i, j] = (rho * (1 / dt * 
                              ((u[i, j+1] - u[i, j-1]) / (2 * dx) + 
                               (v[i+1, j] - v[i-1, j]) / (2 * dy)) -
                              ((u[i, j+1] - u[i, j-1]) / (2 * dx))**2 -
                              2 * ((u[i+1, j] - u[i-1, j]) / (2 * dy) *
                                   (v[i, j+1] - v[i, j-1]) / (2 * dx))-
                              ((v[i+1, j] - v[i-1, j]) / (2 * dy))**2))

    p = pressure_poisson(p, dx, dy, b)

    for i in range(1, Y - 1):
        for j in range(1, X - 1):
            u[i, j] = (u[i, j] -
                       u[i, j] * dt / dx * (u[i, j] - u[i, j-1]) -
                       v[i, j] * dt / dy * (u[i, j] - u[i-1, j]) -
                       dt / (2 * rho * dx) * (p[i, j+1] - p[i, j-1]) +
                       nu * (dt / dx**2 *
                             (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) +
                             dt / dy**2 *
                             (u[i+1, j] - 2 * u[i, j] + u[i-1, j])))
            v[i, j] = (v[i, j] -
                       u[i, j] * dt / dx * (v[i, j] - v[i, j-1]) -
                       v[i, j] * dt / dy * (v[i, j] - v[i-1, j]) -
                       dt / (2 * rho * dy) * (p[i+1, j] - p[i-1, j]) +
                       nu * (dt / dx**2 *
                             (v[i, j+1] - 2 * v[i, j] + v[i, j-1]) +
                             dt / dy**2 *
                             (v[i+1, j] - 2 * v[i, j] + v[i-1, j])))

    # Update boundary conditions
    u[0, :]  = 0
    u[:, 0]  = 0
    v[0, :]  = 0
    v[:, 0]  = 0

    u[-1, :] = 1 
    u[:, -1] = 0
    v[-1, :] = 0
    v[:, -1] = 0

    speed = np.sqrt(u**2 + v**2)
    ax.clear()
    contour = ax.contourf(nX, nY, p, cmap='jet')
    streamplot = ax.streamplot(nX, nY, u, v, color=speed, linewidth=0.5, cmap='jet')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_xlim([0, 2])
    ax.set_ylim([0, 2.25])
    ax.set_title('Pressure and Velocity Field at Time = {}'.format(n))

    return contour, streamplot.lines

ani = animation.FuncAnimation(fig, animate, frames=T, interval=0.5)

plt.show()
