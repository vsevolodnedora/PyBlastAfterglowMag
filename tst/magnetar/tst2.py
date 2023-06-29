import numpy as np
import matplotlib.pyplot as plt

# Constants
gamma = 5/3  # Ratio of specific heats for an ideal gas
delta_t = 0.01  # Time step size

# Function to calculate pressure gradient
def calculate_pressure_gradient(m, v, d, e):
    n_layers = len(m)
    pressure_gradient = np.zeros(n_layers)

    for i in range(1, n_layers):
        delta_v = v[i] - v[i-1]
        delta_e = e[i] - e[i-1]
        delta_m = m[i] - m[i-1]

        pressure_gradient[i] = (delta_v * (delta_e + delta_m * v[i]) +
                                (gamma - 1) * (v[i]**2) * delta_m) / d[i]

    return pressure_gradient

# Function to calculate change in radii
def calculate_radii_change(pressure_gradient, m, v, d):
    n_layers = len(m)
    radii_change = np.zeros(n_layers)

    for i in range(1, n_layers):
        delta_p = pressure_gradient[i] - pressure_gradient[i-1]
        delta_m = m[i] - m[i-1]

        radii_change[i] = -delta_p * d[i] / (2 * delta_m * v[i])

    return radii_change

# Input parameters
m = np.array([1.0, 2.0, 3.0])  # Masses of the layers
v = np.array([0.0, 1.0, 2.0])  # Initial velocities of the layers
d = np.array([0.1, 0.2, 0.3])  # Thicknesses of the layers
e = np.array([10.0, 20.0, 30.0])  # Internal energies of the layers

# Initialize arrays for plotting
n_layers = len(m)
time_steps = 1000
time = np.arange(time_steps) * delta_t
radii = np.zeros((n_layers, time_steps))

# Calculate pressure gradient and radii change over time
for t in range(time_steps):
    pressure_gradient = calculate_pressure_gradient(m, v, d, e)
    radii_change = calculate_radii_change(pressure_gradient, m, v, d)

    # Update radii
    radii[:, t] = radii[:, t-1] + radii_change * delta_t

    # Update velocities and internal energies for next time step
    v += pressure_gradient * delta_t
    e += (gamma - 1) * (v**2) * m * delta_t

# Plotting
plt.figure(figsize=(12, 6))

# Plot pressure gradient
plt.subplot(1, 2, 1)
plt.plot(time, pressure_gradient)
plt.xlabel('Time')
plt.ylabel('Pressure Gradient')
plt.title('Pressure Gradient vs Time')

# Plot change in radii
plt.subplot(1, 2, 2)
for i in range(n_layers):
    plt.plot(time, radii[i])
plt.xlabel('Time')
plt.ylabel('Change in Radii')
plt.title('Change in Radii vs Time')

plt.tight_layout()
plt.show()
