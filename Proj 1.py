import numpy as np
import sys
import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('aspirin.xyz', skiprows=2, usecols=(1, 2, 3))
coord_x = data[:, 0]
coord_y = data[:, 1]
coord_z = data[:, 2]

data_str = np.loadtxt('aspirin.xyz', skiprows=2, dtype=str)
atom = data_str[:, 0]
coord = data_str[:, 1:].astype(float)

#Determine linearity:
#plot all points on a 3d graph
linearity = True #default is true

fig = plot.figure()
axis = fig.add_subplot(111, projection='3d')
axis.scatter(coord_x, coord_y, coord_z)

#take three random atoms and connect a 2D plane so that it intersects all 3 points. 

p1, p2, p3 = data[random.sample(range(len(data)), 3)]

vector1 = np.array(p3) - np.array(p1) # find plane equation using multivariable alg.
vector2 = np.array(p2) - np.array(p1)
cp = np.cross(vector1, vector2)
a, b, c = cp
dp = np.dot(cp, p3)

#stop at the first instance of an atom not on the plane

for coordinates in data:
    if a*coordinates[0] + b*coordinates[1] + c*coordinates[2] != dp:
        linearity = False
        break

print("Linearity: ", linearity)
for atom, coord in zip(atom, coord):
    print(f"Atom: {atom}, Coordinates: {coord}")

plot.show()