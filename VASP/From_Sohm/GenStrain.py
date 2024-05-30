"""
Generates strained POSCAR files.
Generates 6 files corresponding to the 6 strin components.
"""
import numpy as np
import subprocess
import ase.io


SBN = ase.io.read('POSCAR')
pos = SBN.get_positions()
cell = SBN.get_cell()
mag = 1.001



new_cell = np.array(cell).copy()
new_cell[0,0] = new_cell[0,0] * mag
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'xx'])
ase.io.write('xx/POSCAR', SBN, format='vasp')


new_cell = np.array(cell).copy()
new_cell[1,1] = new_cell[1,1] * mag
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'yy'])
ase.io.write('yy/POSCAR', SBN, format='vasp')


new_cell = np.array(cell).copy()
new_cell[2,2] = new_cell[2,2] * mag
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'zz'])
ase.io.write('zz/POSCAR', SBN, format='vasp')


new_cell = np.array(cell).copy()
new_cell[0,1] = new_cell[0,0] * (mag-1)
new_cell[1,0] = new_cell[1,1] * (mag-1)
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'xy'])
ase.io.write('xy/POSCAR', SBN, format='vasp')


new_cell = np.array(cell).copy()
new_cell[0,2] = new_cell[0,0] * (mag-1)
new_cell[2,0] = new_cell[2,2] * (mag-1)
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'xz'])
ase.io.write('xz/POSCAR', SBN, format='vasp')


new_cell = np.array(cell).copy()
new_cell[1,2] = new_cell[1,1] * (mag-1)
new_cell[2,1] = new_cell[2,2] * (mag-1)
SBN.set_cell(new_cell)
subprocess.run(['mkdir', 'yz'])
ase.io.write('yz/POSCAR', SBN, format='vasp')
