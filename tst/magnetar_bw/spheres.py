import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the two points on the inner and outer sphere
phi1, theta1, R1 = 1, 1, 1
phi3, theta3, R3 = 0, 4, 10

# Define the radius of the middle sphere
R2 = (R1 + R3) / 2

# Convert spherical coordinates to Cartesian coordinates
x1 = R1 * np.sin(phi1) * np.cos(theta1)
y1 = R1 * np.sin(phi1) * np.sin(theta1)
z1 = R1 * np.cos(phi1)

x3 = R3 * np.sin(phi3) * np.cos(theta3)
y3 = R3 * np.sin(phi3) * np.sin(theta3)
z3 = R3 * np.cos(phi3)

# Calculate the direction vector of the line between the two points
dx = x3 - x1
dy = y3 - y1
dz = z3 - z1

# Calculate the intersection point of the line with the middle sphere
a = dx**2 + dy**2 + dz**2
b = 2 * (x1*dx + y1*dy + z1*dz)
c = x1**2 + y1**2 + z1**2 - R2**2
disc = b**2 - 4*a*c
t1 = (-b - np.sqrt(disc)) / (2*a)
t2 = (-b + np.sqrt(disc)) / (2*a)
x = x1 + t2*dx
y = y1 + t2*dy
z = z1 + t2*dz

r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arctan2(y, x)
phi = np.arccos(z / r)

print(f"Polar coordinates of intersection point: r={r}, theta={theta}, phi={phi}")


# Plot the three spheres and the intersection point
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
x_inner = R1 * np.sin(v) * np.cos(u)
y_inner = R1 * np.sin(v) * np.sin(u)
z_inner = R1 * np.cos(v)
ax.plot_wireframe(x_inner, y_inner, z_inner, color='r', alpha=0.1)

x_middle = R2 * np.sin(v) * np.cos(u)
y_middle = R2 * np.sin(v) * np.sin(u)
z_middle = R2 * np.cos(v)
ax.plot_wireframe(x_middle, y_middle, z_middle, color='g', alpha=0.1)

x_outer = R3 * np.sin(v) * np.cos(u)
y_outer = R3 * np.sin(v) * np.sin(u)
z_outer = R3 * np.cos(v)
# ax.plot_wireframe(x_outer, y_outer, z_outer, color='b', alpha=0.3)

ax.scatter(x1, y1, z1, color='r', marker='o', label="1")
ax.scatter(x3, y3, z3, color='b', marker='o', label="2")
ax.scatter(x, y, z, color='g', marker='x', label="middle")
ax.plot([x1, x3], [y1, y3], [z1, z3], color='m')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.show()


