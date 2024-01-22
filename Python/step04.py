
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

########################################################
# Step 4: 1D Burgers' Equation
########################################################

X  = 101;                     # Number of spatial points
T  = 100;                     # Number of time steps
v  = 0.07;                    # Diffusion coefficient
dx = 2. * np.pi / (X - 1);    # Spatial step size
dt = dx * v;                  # Time step size

# Initialise spatial grid and initial condition
x  = np.linspace(0, 2*np.pi, num=X)
u  = np.zeros(X)
uA = np.zeros(X)

def AnalyticalSolution(t, v, x):
    phi  = np.exp(-(4*t - x)**2/(4*v*(1 + t))) + np.exp(-(2*np.pi + 4*t - x)**2/(4*v*(1 + t)))
    dphi = (np.exp(-(2*np.pi + 4*t - x)**2/(4*v*(1 + t)))*(4*np.pi + 8*t - 2*x))/(4*v*(1 + t)) + (np.exp(-(4*t - x)**2/(4*v*(1 + t)))*(8*t - 2*x))/(4*v*(1 + t))
    u    = 4 - (2*dphi*v)/phi
    return u

for i in range(X):
    u[i]  = AnalyticalSolution(0, v, x[i])
    uA[i] = AnalyticalSolution(100*dt, v, x[i])

# Setup the plot
fig, ax = plt.subplots()
line, = ax.plot(x, u, 'r', linewidth=1)
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
ax.set_title('1D Burgers\' Equation')

def animate(n):
    global u
    un = u.copy()  
    for i in range(1, X-1):  
        u[i] = (v * dt/dx**2 * (un[i+1] - 2*un[i] + un[i-1])) - (un[i] * dt/dx * (un[i] - un[i-1])) + un[i]
    u[0] = (v * dt / (dx * dx) * (un[1] - 2*un[0] + un[X-2])) - (un[0] * dt / dx * (un[0] - un[X-2])) + un[0]
    u[X-1] = u[0]
    line.set_ydata(u)
    return line,

ani = animation.FuncAnimation(fig, animate, frames=T, interval=50)

plt.show()



