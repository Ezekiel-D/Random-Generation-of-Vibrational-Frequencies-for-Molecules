import numpy as np
import sys
import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D
#_______________________________________________
file = 'Capsaicin.xyz'
linearity = False
#_______________________________________________
prompt_initial = input("use atom and linearity override? (y/n)")    
if prompt_initial == "y":
    N_atom = int(input("Number of atoms in molecule: (int>2)"))
    linearity_override = input("Linear? (y/n)")
    if linearity_override == "y":
        linearity = True
    else:
        linearity = False

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

def determine_vib_df_manual(): 
    if linearity == True:
        vib_df = 3 * N_atom - 5
    else:
        vib_df = 3 * N_atom - 6
    return(vib_df)

def determine_vib_df_file(): 
    if linearity == True:
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
NOTE: There is about a {(N_atom*100)/3.2e+16}% chance of a duplicate item in the Generated vibrational-frequency list.
""")
