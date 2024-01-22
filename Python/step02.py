
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

########################################################
# Step 2: 1D Nonlinear Convection
########################################################

X  = 41;                    # Number of spatial points
T  = 30;                    # Number of time steps
dx = 2. / (X - 1.);         # Spatial step size
dt = 0.2 * dx;              # Time step size

# Initialise spatial grid and initial condition
x = np.linspace(0, 5, num=X)
u = np.ones(X)
u[(x >= 0.5) & (x <= 1)] = 2

# Setup the plot
fig, ax = plt.subplots()
line, = ax.plot(x, u, 'r', linewidth=1)
ax.set_xlim(0, 2)
ax.set_ylim(0.5, 2.5)
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
ax.set_title('1D Nonlinear Convection')

def animate(n):
    global u
    un = u.copy()  
    for i in range(1, X):  
        u[i] = un[i] - un[i] * (un[i] - un[i-1]) * dt / dx
    line.set_ydata(u)
    return line,

ani = animation.FuncAnimation(fig, animate, frames=T, interval=50)

plt.show()



