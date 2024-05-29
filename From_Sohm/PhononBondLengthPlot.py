"""
Plots bond length of selected pair of atoms as a function of phonon eigendisplacements.
"""
import numpy as np
import ase.io
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
# import sys
# sys.path.append('../')
from AnaddbParser import AnaddbParser



anaddb = AnaddbParser("../SBN_anaddb.abo")
EO_tensors = anaddb.get_EO()
r_33 = [i[2,2] for i in EO_tensors["EO Tensor"]]


natoms = 45
amplitude = 1
pair = (7, 37)


sbn = ase.io.read("POSCAR")
cell = np.diag(np.array(sbn.cell))
pos = sbn.get_positions()
dist_orig = np.abs(pos[pair[0]-1] - pos[pair[1]-1])
for i in range(3):
    if dist_orig[i] > cell[i]/2:
        dist_orig[i] = cell[i]-dist_orig[i]
dist_orig = np.linalg.norm(dist_orig)

df = pd.DataFrame(columns=["Mode", "Change in bond length"])


for mode in np.arange(4,3*natoms+1):
    with open("modulation.conf", 'w') as file:
        file.write(f"FORCE_CONSTANTS = READ\nDIM = 1 1 3\nMODULATION = 1 1 1, 0.0 0.0 0.0 {mode} {amplitude} 0")

    subprocess.run("phonopy -q modulation.conf", shell=True, executable='/bin/bash')

    sbn = ase.io.read("MPOSCAR")
    pos = sbn.get_positions()
    dist = np.abs(pos[pair[0]-1] - pos[pair[1]-1])
    for i in range(3):
        if dist[i] > cell[i]/2:
            dist[i] = cell[i]-dist[i]
    dist = np.linalg.norm(dist)
    delta_dist = np.abs(dist_orig - dist)

    df.loc[len(df.index)] = [mode, delta_dist]


df = df.set_index("Mode")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
fig.set_figwidth(15)
fig.set_figheight(5)

ax1.plot(EO_tensors["Mode"], r_33, 'b-')
ax2.plot(EO_tensors["Mode"], df["Change in bond length"],'r-')
ax1.set_xlabel("phonon mode number")
ax1.set_ylabel("r_33 (pm/V)", color="b")
ax2.set_ylabel("relative Nb-O bond length", color="r")
ax1.set_title("Mode Decomposed EO Tensor")
plt.savefig("BondLengthPlot.pdf")

