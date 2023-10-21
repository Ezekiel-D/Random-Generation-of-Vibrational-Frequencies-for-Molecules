import numpy as np
import sys
import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D
#_______________________________________________
file = 'benzene.xyz'
#_______________________________________________
linearity = input(f"State linearity of {file}: (t/f)")
if linearity == "t":
    linearity = True
else:
    linearity = False
prompt_initial = input("use atom and linearity override? (y/n)")    
if prompt_initial == "y":
    N_atom = int(input("Number of atoms in molecule: (int>2)"))
    linearity_override = input("Linear? (y/n)")
    if linearity_override == "y":
        linearity = True
    else:
        linearity = False
##Step 5

##Step 1
#extracts x, y, z coordinates
data = np.loadtxt(file, skiprows=2, usecols=(1, 2, 3))
coord = data[:, 0:]
coord_x = data[:, 0]
coord_y = data[:, 1]
coord_z = data[:, 2]

#extracts atom
data_str = np.loadtxt(file, skiprows=2, dtype=str)
atom = data_str[:, 0]


#vibrational degrees freedom acquired from data here
if prompt_initial != "y":
    N_atom = len(atom)
df = 3 * N_atom

##EXT: Determine linearity
#plot all points on a 3d graph

pltfig = plot.figure()
graph = pltfig.add_subplot(111, projection='3d')
graph.scatter(coord_x, coord_y, coord_z)

planarity = "N/A" #default
def planar_test():
#take three random atoms and connect a 2D plane so that it intersects all 3 points. 
    p1, p2, p3 = data[random.sample(range(len(data)), 3)]
#find plane equation using multivariable alg.
    vector1 = np.array(p3) - np.array(p1) 
    vector2 = np.array(p2) - np.array(p1)
    cp = np.cross(vector1, vector2)
    a, b, c = cp
    d = np.dot(cp, p3)
#plot plane
    x = np.linspace(min(data[:, 0]), max(data[:, 0]))
    y = np.linspace(min(data[:, 1]), max(data[:, 1]))
    X, Y = np.meshgrid(x, y)
    Z = (d - a*X - b*Y) / c #ax + by + cz = d
    graph.plot_surface(X, Y, Z, alpha=0.1)
#stop at the first instance of an atom not on the plane
    for coordinates in data:
        if not np.isclose(a*coordinates[0] + b*coordinates[1] + c*coordinates[2], d, rtol=1e-02, atol=1e-02):
            planarity = False
            return(planarity)
    planarity = True
    return(planarity)
     

if not linearity:
    planarity = planar_test()

def determine_vib_df_manual(): 
    if linearity:
        vib_df = 3 * N_atom - 5
    else:
        vib_df = 3 * N_atom - 6
    return(vib_df)

def determine_vib_df_file(): 
    if linearity:
        vib_df = df - 5
    else:
        vib_df = df - 6
    return(vib_df)

if prompt_initial == "y":
    vib_df = determine_vib_df_manual()
else:
    vib_df = determine_vib_df_file()


##Step 2
#generating vibrational_list
vibrational_list = [] 
for i in range(vib_df):
    vibrational_list.append(np.random.uniform(0,3200))
vibrational_list.sort()

#sorting and creating subarrays from vibrational_spectrum
lowfreq = vibrational_list[:] #in the case that all generated frequencies are below 800
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
N_given = np.loadtxt(file, skiprows=0, usecols=(0), max_rows=1)  

if prompt_initial == "y":
    Nmatch = "NA"
elif N_given == N_atom:
    Nmatch = True

##Step 4
print(f"""
    Number of atoms: {N_atom} (Match: {Nmatch})
    Linearity: {linearity}
    Planarity: {planarity}
    Vibrational degrees of freedom: {vib_df}
""")
if prompt_initial != "y":
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
NOTE: There is about a {(N_atom*100)/3.2e+16}% chance of a duplicate item in the generated vibrational-frequency list.
""")

if prompt_initial != "y":
    prompt = input("View plot of atoms from coordinate file? (y/n)")
    if prompt == "y":
        plot.show()
    else:
        StopIteration
