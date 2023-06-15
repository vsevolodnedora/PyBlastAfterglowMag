import numpy as np
import matplotlib.pyplot as plt

# Circle parameters
r1 = 2  # Radius of inner circle
r2 = 4  # Radius of middle circle
r3 = 6  # Radius of outer circle

# Point coordinates
theta1 = np.pi / 6  # Angle for point on inner circle (in radians)
theta3 = np.pi / 3  # Angle for point on outer circle (in radians)

# Convert polar coordinates to Cartesian coordinates
x1, y1 = r1 * np.cos(theta1), r1 * np.sin(theta1)
x3, y3 = r3 * np.cos(theta3), r3 * np.sin(theta3)

# Parametric equation of the line
x = np.linspace(x1, x3, num=100)
y = np.linspace(y1, y3, num=100)

# Find the intersection point with the middle circle
t = (r2 - r1) / (r3 - r1)  # Parameter for the intersection point
x_intersection = x1 + t * (x3 - x1)
y_intersection = y1 + t * (y3 - y1)

# Plotting
circle1 = plt.Circle((0, 0), r1, color='r', fill=False)
circle2 = plt.Circle((0, 0), r2, color='g', fill=False)
circle3 = plt.Circle((0, 0), r3, color='b', fill=False)
plt.gca().add_patch(circle1)
plt.gca().add_patch(circle2)
plt.gca().add_patch(circle3)
plt.plot(x, y, 'k--', label='Line')
plt.plot(x_intersection, y_intersection, 'ro', label='Intersection')
plt.plot([x1, x3], [y1, y3], 'ko', label='Given Points')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
plt.legend()
plt.title('Intersection of Line with Middle Circle')
plt.grid(True)
plt.show()