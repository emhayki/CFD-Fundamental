
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

########################################################
# Step 12: Channel Flow
########################################################

X = 21
Y = 21
T = 50
F = 1

rho = 1
dt  = 0.01
nu  = 0.10

xmin = 0
xmax = 2

ymin = 0
ymax = 1

dx = (xmax - xmin) / (X - 1)
dy = (ymax - ymin) / (Y - 1)

p, b = np.zeros((Y, X)), np.zeros((Y, X))
u, v = np.zeros((Y, X)), np.zeros((Y, X))

x = np.linspace(0, 2, X)
y = np.linspace(0, 2, Y)
nX, nY = np.meshgrid(x, y)

# Lid-driven cavity condition
u[-1, :] = 1 

# Setup the plot
fig, ax = plt.subplots()
contour = ax.contourf(nX, nY,  p, cmap='jet')
quiver  = ax.quiver(nX, nY, u, v, cmap='jet')

# Pressure Poisson equation
def pressure_poisson(p, dx, dy, b):
    for _ in range(50):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) / (2 * (dx**2 + dy**2)) - dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

        # Boundary Condition
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 + (pn[2:, -1] + pn[0:-2, -1]) * dx**2) / (2 * (dx**2 + dy**2)) - dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])

        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 + (pn[2:, 0] + pn[0:-2, 0]) * dx**2) / (2 * (dx**2 + dy**2)) - dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])
        
        # Wall Boundary Condition
        p[-1, :] = p[-2, :]  
        p[ 0, :] = p[ 1, :]  
    return p


# Animation function
def animate(n):

    global u, v, p

    un = u.copy()
    vn = v.copy()

    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    
    # Boundary Condition
    b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx) + (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) - ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 - 2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) * (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) - ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))
    b[1:-1,  0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) + (v[2:,  0] - v[0:-2,  0]) / (2 * dy)) - ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 - 2 * ((u[2:,  0] - u[0:-2,  0]) / (2 * dy) * (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx)) - ((v[2:,  0] - v[0:-2,  0]) / (2 * dy))**2))
    
    p = pressure_poisson(p, dx, dy, b)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) - vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1]) - dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu * (dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) +  F * dt)
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) - vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu * (dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) + dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # Boundary Condition
    u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * (un[1:-1, -1] - un[1:-1, -2]) - vn[1:-1, -1] * dt / dy * (un[1:-1, -1] - un[0:-2, -1]) - dt / (2 * rho * dx) * (p[1:-1, 0] - p[1:-1, -2]) + nu * (dt / dx**2 * (un[1:-1, 0] - 2 * un[1:-1, -1] + un[1:-1, -2]) + dt / dy**2 * (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)
    v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx * (vn[1:-1, -1] - vn[1:-1, -2]) - vn[1:-1, -1] * dt / dy * (vn[1:-1, -1] - vn[0:-2, -1]) - dt / (2 * rho * dy) * (p[2:,  -1] - p[0:-2, -1]) + nu * (dt / dx**2 * (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) + dt / dy**2 * (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx * (un[1:-1, 0] - un[1:-1, -1]) - vn[1:-1, 0] * dt / dy * (un[1:-1, 0] - un[0:-2, 0]) - dt / (2 * rho * dx) * (p[1:-1, 1] - p[1:-1, -1]) + nu * (dt / dx**2 * (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) + dt / dy**2 * (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)
    v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx * (vn[1:-1, 0] - vn[1:-1, -1]) - vn[1:-1, 0] * dt / dy * (vn[1:-1, 0] - vn[0:-2, 0]) - dt / (2 * rho * dy) * (p[2:,   0] - p[0:-2,  0]) + nu * (dt / dx**2 * (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) + dt / dy**2 * (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))

    # Wall Boundary Condition
    u[0, :] = 0
    v[0, :] = 0

    u[-1, :] = 0
    v[-1, :] = 0

    ax.clear()
    contour = ax.contourf(nX, nY, p, cmap='jet')
    quiver = ax.quiver(nX, nY, u, v, cmap='jet')
    ax.set_xlim([0, 2])
    ax.set_ylim([0, 2])
    return contour, quiver

ani = animation.FuncAnimation(fig, animate, frames=T, interval=0.5)

plt.show()
