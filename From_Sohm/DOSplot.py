"""
Plots DOS as calculated by VASP.
"""
import pyprocar

pyprocar.dosplot(mode='stack_species', savefig='dos.pdf', title='Total Density of States')
