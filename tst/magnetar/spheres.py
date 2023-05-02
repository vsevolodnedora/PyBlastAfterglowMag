import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the two points on the inner and outer sphere
phi1, theta1, R1 = 1, 1, 1
phi2, theta2, R3 = 0, 4, 10

# Define the radius of the middle sphere
R2 = (R1 + R3) / 2

# Convert spherical coordinates to Cartesian coordinates
x1 = R1 * np.sin(phi1) * np.cos(theta1)
y1 = R1 * np.sin(phi1) * np.sin(theta1)
z1 = R1 * np.cos(phi1)

x2 = R3 * np.sin(phi2) * np.cos(theta2)
y2 = R3 * np.sin(phi2) * np.sin(theta2)
z2 = R3 * np.cos(phi2)

# Calculate the direction vector of the line between the two points
dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

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

ax.scatter(x1, y1, z1, color='r', marker='o')
ax.scatter(x2, y2, z2, color='b', marker='o')
ax.scatter(x, y, z, color='g', marker='x')
ax.plot([x1, x2], [y1, y2], [z1, z2], color='m')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()





exit(0)

import numpy as np
import matplotlib.pyplot as plt

# Define the parameters of the problem
R1 = 1.0    # Radius of inner sphere
R2 = 2.0    # Radius of outer sphere
theta1 = np.pi/4    # Angle theta of point on inner sphere
phi1 = np.pi/3      # Angle phi of point on inner sphere
theta3 = np.pi/2    # Angle theta of point on outer sphere
phi3 = 0.0          # Angle phi of point on outer sphere

# Compute the coordinates of the two points
x1 = R1 * np.sin(theta1) * np.cos(phi1)
y1 = R1 * np.sin(theta1) * np.sin(phi1)
z1 = R1 * np.cos(theta1)

x3 = R2 * np.sin(theta3) * np.cos(phi3)
y3 = R2 * np.sin(theta3) * np.sin(phi3)
z3 = R2 * np.cos(theta3)

# Compute the coordinates of the intersection point
t = (R2 - R1) / np.sqrt((x3 - x1)**2 + (y3 - y1)**2 + (z3 - z1)**2)
xi = x1 + t*(x3 - x1)
yi = y1 + t*(y3 - y1)
zi = z1 + t*(z3 - z1)

# Plot the spheres and the points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xs1 = R1*np.sin(v)*np.cos(u)
ys1 = R1*np.sin(v)*np.sin(u)
zs1 = R1*np.cos(v)
ax.plot_wireframe(xs1, ys1, zs1, color="b")

xs2 = R2*np.sin(v)*np.cos(u)
ys2 = R2*np.sin(v)*np.sin(u)
zs2 = R2*np.cos(v)
ax.plot_wireframe(xs2, ys2, zs2, color="r")

ax.scatter([x1, xi], [y1, yi], [z1, zi], color="g")
ax.scatter([x3], [y3], [z3], color="k")

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.plot([x1, x2], [y1, y2], [z1, z2], color='m')
plt.show()
