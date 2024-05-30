"""
Plots band structure (or band+DOS) calculated by VASP.
"""
import pyprocar
import numpy as np

# atoms = np.arange(16,46)
# pyprocar.bandsplot('PROCAR',outcar='OUTCAR',elimit=[-10,10],mode='parametric',atoms=atoms,code='vasp',savefig='bands.pdf',kpointsfile='KPOINTS')

# pyprocar.bandsdosplot(bands_file='PROCAR',dos_file='vasprun.xml',outcar='OUTCAR',elimit=[-10,10],dos_spins=[0],code='vasp',savefig='bandsdos.pdf',kpointsfile='KPOINTS')

pyprocar.bandsplot('PROCAR',outcar='OUTCAR',elimit=[-10,10],mode='parametric',code='vasp',savefig='bands.pdf',kpointsfile='KPOINTS')
