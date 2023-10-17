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
coord = data_str[:, 1].astype(float)

N_atom = len(atom)
#Determine linearity:
#plot all points on a 3d graph
linearity = True #default is true

fig = plot.figure()
graph = fig.add_subplot(111, projection='3d')
graph.scatter(coord_x, coord_y, coord_z)

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

df = 3 * N_atom #degrees freedom
vib_df = df - 6


vibrational_spectrum = [] #generating vibrational_spectrum
for i in range(vib_df + 1):
    vibrational_spectrum.append(np.random.uniform(0,3200))
vibrational_spectrum.sort()
print(vibrational_spectrum)
for i in range(vib_df + 1):
    if (vibrational_spectrum[i] > 800):
        lowfreq = vibrational_spectrum[:i]
        pointer = i
        break
for i in range(pointer, vib_df):
    if (vibrational_spectrum[i] > 1600):
        fingerprint = vibrational_spectrum[pointer:i]
        highfreq = vibrational_spectrum[i:]
        break

N_lo = len(lowfreq)
N_hi = len(highfreq)

print(f"Number of atoms: {N_atom}")  
print(f"Vibrational df: {vib_df}")
print(f"Linearity: {linearity}")
for atom, coord in zip(atom, coord):
    print(f"Atom: {atom}, Coordinates: {coord}")
print(f"Low Frequency domain ({N_lo} atoms): {lowfreq}")
print(f"Fingerprint domain ({N_hi-N_lo} atoms): {fingerprint}")
print(f"High Frequency domain ({N_hi} atoms): {highfreq}")

