import numpy as np
import matplotlib.pyplot as plt

# Parameters
E0 = 1e53  # Initial energy (erg)
Gamma0 = 1000  # Initial Lorentz factor
theta0 = 0.2  # Opening angle (radians)
n = 3  # Power-law index for the density profile

# Constants
c = 3e10  # Speed of light (cm/s)
Msun = 1.989e33  # Solar mass (g)
pc = 3.086e18  # Parsec (cm)

# Conversion factors
erg_to_eV = 6.242e11  # erg to eV conversion
cm_to_pc = 3.241e-19  # cm to pc conversion

# Function to calculate density profile
def density_profile(r, t):
    t_dyn = (theta0 * c) / (Gamma0 * c**2)
    r_dyn = Gamma0**2 * c**2 * t_dyn
    eta = r / r_dyn
    rho = (3 * E0 * n * (n - 3) / (4 * np.pi * (n - 2)**2 * r_dyn**3)) * (1 + t / t_dyn)**(-3 * n / (n - 2))
    rho *= (1 + eta**2)**(-(n + 2) / 2)
    return rho

# Time steps for plotting
timesteps = np.logspace(0, 4, num=5)  # Modify the range and number of time steps as needed

# Distance range for plotting
r_min = 0.01 * c * theta0 / (Gamma0**2 * c**2)
r_max = 10 * c * theta0 / (Gamma0**2 * c**2)
r = np.logspace(np.log10(r_min), np.log10(r_max), num=1000)

# Plotting
for t in timesteps:
    rho = density_profile(r, t)
    plt.loglog(r / (pc * cm_to_pc), rho * erg_to_eV / (Msun * pc**-3), label='t = {:.2e} s'.format(t))

plt.xlabel('Distance (pc)')
plt.ylabel('Density (Msun/pc^3)')
plt.title('Density Profile Behind Relativistic Blast Wave')
plt.legend()
plt.grid(True)
plt.show()
