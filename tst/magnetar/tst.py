import numpy as np
import matplotlib.pyplot as plt

# Constants
dt = 0.01  # Time step
dx = 0.01  # Spatial step
t_max = 10.0  # Maximum simulation time
x_max = 10.0  # Maximum spatial extent

# Grid parameters
Nt = int(t_max / dt) + 1  # Number of time steps
Nx = int(x_max / dx) + 1  # Number of spatial steps

# Arrays to store PWN and ejecta data
pwn_density = np.zeros(Nx)
pwn_velocity = np.zeros(Nx)
ejecta_density = np.zeros(Nx)
ejecta_velocity = np.zeros(Nx)
x_values = np.linspace(0.0, x_max, Nx)

# Initial conditions
pwn_density[:] = 1.0
ejecta_density[:] = 1.0

# Main simulation loop
for t in range(Nt):
    # Update PWN density and velocity
    pwn_density[1:-1] += dt * (-pwn_velocity[1:-1] * np.diff(pwn_density) / dx)
    pwn_velocity[1:-1] += dt * (-pwn_velocity[1:-1] * np.diff(pwn_velocity) / dx)

    # Update ejecta density and velocity
    ejecta_density[1:-1] += dt * (-ejecta_velocity[1:-1] * np.diff(ejecta_density) / dx)
    ejecta_velocity[1:-1] += dt * (-ejecta_velocity[1:-1] * np.diff(ejecta_velocity) / dx)

    # Reflective boundary conditions
    pwn_density[0] = pwn_density[1]
    pwn_density[-1] = pwn_density[-2]
    pwn_velocity[0] = -pwn_velocity[1]
    pwn_velocity[-1] = -pwn_velocity[-2]

    ejecta_density[0] = ejecta_density[1]
    ejecta_density[-1] = ejecta_density[-2]
    ejecta_velocity[0] = -ejecta_velocity[1]
    ejecta_velocity[-1] = -ejecta_velocity[-2]

# Plotting the results
plt.figure(figsize=(10, 5))
plt.plot(x_values, pwn_density, label='PWN Density')
plt.plot(x_values, ejecta_density, label='Ejecta Density')
plt.xlabel('Position')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
plt.title('Evolution of PWN and Ejecta Shells')
plt.show()
