"""
Generates displacements of ions along a phonon eigenvector. 
Displacement magnitudes are entered as command line arguments.
"""
import numpy as np
import yaml
import sys
import subprocess
import ase.io


# Get phonon frequencies and eigenvectors from Phonopy.
def get_phonons(qpoint_file = 'qpoints.yaml', band=1):
    with open(qpoint_file, 'r') as file:
        data = yaml.safe_load(file)

    freq = data['phonon'][0]['band'][band-1]['frequency']     # frequency of band
    eig = np.array(data['phonon'][0]['band'][band-1]['eigenvector'])      # eigenvector of band
    eig = eig[:,:,0]
    return (freq, eig)

freq, eig = get_phonons(band=4)

# Loading structure data with ASE.
SBN = ase.io.read('POSCAR')
pos = SBN.get_positions()
cell = SBN.get_cell()
mass = SBN.get_masses()

eig = eig / np.sqrt(np.tile(mass, (3,1)).T)     # eigen-displacements = eigenvectors / sqrt(mass).

# Reading displacement magnitudes.
mag = sys.argv[1:]
mag = [float(i) for i in mag]

# Writing displaced structures in POSCAR files.
for i in mag:
    pos_new = pos + i * eig
    pos_new[:,0] = np.mod(pos_new[:,0], cell[0,0])
    pos_new[:,1] = np.mod(pos_new[:,1], cell[1,1])
    pos_new[:,2] = np.mod(pos_new[:,2], cell[2,2])

    SBN.set_positions(pos_new)
    if i >= 0:
        subprocess.run(['mkdir', str(i)])
        ase.io.write(f'{i}/POSCAR', SBN, format='vasp')
    else:
        subprocess.run(['mkdir', 'm'+str(-i)])
        ase.io.write(f'm{-i}/POSCAR', SBN, format='vasp')
    
