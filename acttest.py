import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

# Load atomic coordinates from an .xyz file
data = np.loadtxt('Benzene.xyz', skiprows=2, usecols=(1, 2, 3))

# Plot all points on a 3D graph
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data[:, 0], data[:, 1], data[:, 2])

# Select three random atoms
indices = random.sample(range(data.shape[0]), 3)
points = data[indices, :]

#if np.linalg.matrix_rank(points) < 3:
   #points = data[indices, :]
# Calculate the equation of the plane
p1, p2, p3 = points
v1 = p3 - p1
v2 = p2 - p1
cp = np.cross(v1, v2)
a, b, c = cp
d = np.dot(cp, p3)

# Generate x and y coordinates for the plane
x = np.linspace(min(data[:, 0]), max(data[:, 0]), 10)
y = np.linspace(min(data[:, 1]), max(data[:, 1]), 10)
X, Y = np.meshgrid(x, y)

# Calculate z coordinates for the plane
Z = (d - a*X - b*Y) / c

# Plot the plane
ax.plot_surface(X, Y, Z, alpha=0.5, rstride=100, cstride=100)

planarity = True
print(planarity)
# Check if other atoms lie on the plane
for coord in data:
    if not np.isclose(a*coord[0] + b*coord[1] + c*coord[2], d, rtol=1e-02, atol=1e-02):
        planarity = False
        break
print(planarity)
