"""
Generate POSCAR files with atomic displacements.
"""
import numpy as np
import yaml
import sys
import subprocess
import ase.io


SBN = ase.io.read('POSCAR')
pos = SBN.get_positions()
cell = SBN.get_cell()


atoms = np.arange(45)
directions = [0, 1, 2]
magnitudes = [-0.01, 0.01]

for atom in atoms:
    for dir in directions:
        subprocess.run(['mkdir', f'{atom*3+dir}'])

        for mag in magnitudes:
            new_pos = pos.copy()
            new_pos[atom][dir] += mag
            new_pos[:,0] = np.mod(new_pos[:,0], cell[0,0])
            new_pos[:,1] = np.mod(new_pos[:,1], cell[1,1])
            new_pos[:,2] = np.mod(new_pos[:,2], cell[2,2])

            SBN.set_positions(new_pos)

            if mag >= 0:
                subprocess.run(['mkdir', f'{atom*3+dir}/plus'])
                ase.io.write(f'{atom*3+dir}/plus/POSCAR', SBN, format='vasp')
            else:
                subprocess.run(['mkdir', f'{atom*3+dir}/minus'])
                ase.io.write(f'{atom*3+dir}/minus/POSCAR', SBN, format='vasp')
