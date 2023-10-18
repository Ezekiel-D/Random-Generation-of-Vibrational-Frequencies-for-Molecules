import numpy as np
import sys
import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D

##Step 1
#extracts x, y, z coordinates
data = np.loadtxt('aspirin.xyz', skiprows=2, usecols=(1, 2, 3))
coord_x = data[:, 0]
coord_y = data[:, 1]
coord_z = data[:, 2]

#extracts atom and full coordinates
data_str = np.loadtxt('aspirin.xyz', skiprows=2, dtype=str)
atom = data_str[:, 0]
coord = data_str[:, 1:].astype(float)

#Determine linearity:
#plot all points on a 3d graph
linearity = True #default is true

fig = plot.figure()
graph = fig.add_subplot(111, projection='3d')
graph.scatter(coord_x, coord_y, coord_z)

#take three random atoms and connect a 2D plane so that it intersects all 3 points. 
p1, p2, p3 = data[random.sample(range(len(data)), 3)]

#find plane equation using multivariable alg.
vector1 = np.array(p3) - np.array(p1) 
vector2 = np.array(p2) - np.array(p1)
cp = np.cross(vector1, vector2)
a, b, c = cp
dp = np.dot(cp, p3)

#stop at the first instance of an atom not on the plane
for coordinates in data:
    if a*coordinates[0] + b*coordinates[1] + c*coordinates[2] != dp:
        linearity = False
        break

##Step 4
#manually type in number of atoms and linearity state here
N_atoms_override = 20
linearity_overide = False
def determine_vib_df_manual():  #entire number of atoms then True or False
    if linearity_overide == True:
        vib_df = 3 * N_atoms_override - 5
    else:
        vib_df = 3 * N_atoms_override - 6
    return(vib_df)
vib_df_manual = determine_vib_df_manual()
#vibrational degrees freedom acquired from data here
N_atom = len(atom)
df = 3 * N_atom

def determine_vib_df(): 
    if linearity == True:
        vib_df = df - 5
    else:
        vib_df = df - 6
    return(vib_df)

vib_df = determine_vib_df()

##Step 2
#generating vibrational_list
vibrational_list = [] 
for i in range(vib_df):
    vibrational_list.append(np.random.uniform(0,3200))
vibrational_list.sort()

#sorting and creating subarrays from vibrational_spectrum
for i in range(vib_df + 1):
    if (vibrational_list[i] > 800):
        lowfreq = vibrational_list[:i]
        pointer = i
        break
for i in range(pointer, vib_df):
    if (vibrational_list[i] > 1600):
        fingerprint = vibrational_list[pointer:i]
        highfreq = vibrational_list[i:]
        break

##Step 3
#extracts number of atom given by document and sees whether it matches with the 'atom' array
Nmatch = False
N_given = np.loadtxt('aspirin.xyz', skiprows=0, usecols=(0), max_rows=1)  
if (N_given == N_atom):
    Nmatch = True

##Step 4


print(f"""
    Number of atoms: {N_atom} (Match: {Nmatch})
    Linearity: {linearity}
    Vibrational degrees of freedom: {vib_df}
""")

print("Transcribed coordinate file:")
for atom, coord in zip(atom, coord):
    print(f"Atom: {atom}, Coordinates: {coord}")

print(f"""
    Generated vibrational-frequency list ({len(vibrational_list)} items):
{vibrational_list}\n
    Low-Frequency range ({len(lowfreq)} items):
{lowfreq}\n
    Fingerprint range ({len(fingerprint)} items):
{fingerprint}\n
    High-Frequency range ({len(highfreq)} items):
{highfreq}\n
NOTE: There is a less than 3E-12% chance of a duplicate item in the Generated vibrational-frequency list for a 100-atoms molecule.
""")

prompt = input("View plot of atoms from coordinate file? (y/n)")
if prompt == "y":
    plot.show()
else:
    StopIteration
