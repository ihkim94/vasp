"""
Calculates elasto-optic tensor from a finite difference VASP calculation.
Can also calculate piezo Pockels tensor if the piezoelectric tensor is given.
"""
import numpy as np
import pymatgen.io.vasp
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)


piezo = np.array([[0.00000008, -0.00000004, -13.64932563],
                    [-0.00000008, -0.00000002, -12.55065303],
                    [0.00000001, 0.00000008, 45.30773324],
                    [-1.26880611, 11.27689195, 0.00000025],
                    [9.37229998, -4.39903353, 0.00000079],
                    [-0.00000003, 0.00000012, -10.20358875]])



mag = 0.001
elasto_optic = np.zeros((3,3,6))

vasprun = pymatgen.io.vasp.outputs.Vasprun(
            f"../03-Dielectric/cart/vasprun.xml",
            parse_dos=False,
            parse_eigen=False,
            parse_projected_eigen=False)
eps_0 = np.array(vasprun.epsilon_static)


for i in range(6):
    vasprun = pymatgen.io.vasp.outputs.Vasprun(
                f"{i+1}/vasprun.xml",
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False)
    eps = np.array(vasprun.epsilon_static)
    eps = (eps + eps.T) / 2         # Symmetrize epsilon.

    elasto_optic[:,:,i] = (eps - eps_0) / mag


elasto_voigt = np.zeros((6,6))
elasto_voigt[0,:] = elasto_optic[0,0,:]
elasto_voigt[1,:] = elasto_optic[1,1,:]
elasto_voigt[2,:] = elasto_optic[2,2,:]
elasto_voigt[3,:] = (elasto_optic[1,2,:]+elasto_optic[2,1,:])/2
elasto_voigt[4,:] = (elasto_optic[0,2,:]+elasto_optic[2,0,:])/2
elasto_voigt[5,:] = (elasto_optic[0,1,:]+elasto_optic[1,0,:])/2
print(elasto_voigt)



piezo_pockels = np.zeros((3,3,3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            piezo_pockels[i,j,k] =  elasto_optic[i,j,0]*piezo[0,k] + elasto_optic[i,j,1]*piezo[1,k] + elasto_optic[i,j,2]*piezo[2,k] + 2*elasto_optic[i,j,3]*piezo[3,k] + 2*elasto_optic[i,j,4]*piezo[4,k] + 2*elasto_optic[i,j,5]*piezo[5,k]
pockels_voigt = np.zeros((6,3))
pockels_voigt[0,:] = piezo_pockels[0,0,:]
pockels_voigt[1,:] = piezo_pockels[1,1,:]
pockels_voigt[2,:] = piezo_pockels[2,2,:]
pockels_voigt[3,:] = piezo_pockels[1,2,:]
pockels_voigt[4,:] = piezo_pockels[0,2,:]
pockels_voigt[5,:] = piezo_pockels[0,1,:]
print(pockels_voigt)


