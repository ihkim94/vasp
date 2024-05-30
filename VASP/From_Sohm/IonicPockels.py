"""
Calculate Raman susceptibility, polarity, and ionic Pockels tensor.
Based on finite difference calculation in VASP.
"""
import numpy as np
import pymatgen.io.vasp
import xml.etree.ElementTree as ET
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import ase.io
np.set_printoptions(suppress=True)


# Get Born effective charges from vasprun.xml
def get_Born(vasprun_born='vasprun.xml'):
    tree = ET.parse(vasprun_born)
    root = tree.getroot()

    natom = 45
    born = np.zeros((natom,3,3))

    x = root.findall("./calculation/array[@name='born_charges']/set/v")
    for i,j in enumerate(x):
        born[i//3][i%3] = j.text.split()

    return born

# Get phonon frequencies and eigen-displacements from Phonopy
def get_phonons(qpoint_file='qpoints.yaml', band=4, poscar='POSCAR'):
    with open(qpoint_file, 'r') as file:
        data = yaml.safe_load(file)

    freq = data['phonon'][0]['band'][band-1]['frequency']     # frequency of band
    eig = np.array(data['phonon'][0]['band'][band-1]['eigenvector'])      # eigenvector of band
    eig = eig[:,:,0]

    structure = ase.io.read(poscar)
    mass = structure.get_masses()
    eig = eig / np.sqrt(np.tile(mass, (3,1)).T)     # eigendisplacements = eigenvectors / sqrt(mass).

    return (freq, eig)


def calculate_polarity(born, phonon):
    return np.einsum('ijk,ik->j', born, phonon)


def calculate_index(dir, x=0.60, lam=1550):
    """
    Calulate refractive index of SBN.
    Arguments:
        dir: 'o' for ordinary, 'e' for extraordinary
        x: Sr concentration
        lam: wavelength in nm

    Fitting parameters calculated by:
    C. David, A. Tunyagi, K. Betzler, and M. Wöhlecke, Phys. Status Solidi B 244, 2127 (2007).
    """

    if dir == 'o':
        A0 = 5.002; A1 = 0
        B0 = 1.272E+5; B1 = 0
        C0 = 6.261E+4; C1 = 0
        D0 = 5.207E-8; D1 = 0

    if dir == 'e':
        A0 = 4.712; A1 = 0.291
        B0 = 9.242E+4; B1 = 4.398E+4
        C0 = 5.679E+4; C1 = 3.62E+3
        D0 = 6.21E-8; D1 = -2.27E-8

    return np.sqrt(A0 + A1*x + (B0+B1*x)/(lam**2-C0-C1*x) - (D0+D1*x)*lam**2)


def susceptibility_derivatives(disp=0.01, natoms=45):
    d_chi_d_tau = []

    for i in range(natoms*3):
        out_m = pymatgen.io.vasp.outputs.Vasprun(
            f"{i}/minus/vasprun.xml",
            parse_dos=False,
            parse_eigen=False,
            parse_projected_eigen=False,
            parse_potcar_file=False)
        chi_m = (np.array(out_m.epsilon_static) - 1) / (4*np.pi)
        #chi_m = np.array(out_m.epsilon_static)

        out_p = pymatgen.io.vasp.outputs.Vasprun(
            f"{i}/plus/vasprun.xml",
            parse_dos=False,
            parse_eigen=False,
            parse_projected_eigen=False,
            parse_potcar_file=False)
        chi_p = (np.array(out_p.epsilon_static) - 1) / (4*np.pi)
        #chi_p = np.array(out_p.epsilon_static)

        d_chi_d_tau.append((chi_p - chi_m)/(2*disp))

    d_chi_d_tau = np.array(d_chi_d_tau)
    d_chi_d_tau = np.reshape(d_chi_d_tau, (natoms, 3, 3, 3))
    err = np.sum(d_chi_d_tau, axis=0) / natoms
    for i in range(natoms):
        d_chi_d_tau[i] = d_chi_d_tau[i] - err

    return d_chi_d_tau


def calculate_Raman(d_chi_d_tau, phonon):
    return np.einsum('ijkl,ij->kl', d_chi_d_tau, phonon)


def calculate_Pockels_ionic(alphas, polarities, freqs): 
    # ne = calculate_index('e', 0.6, 1550)       # Using wavelength = 1550 nm.
    # no = calculate_index('o', 0.6, 1550)
    # print(ne, no)
    # n = np.array([1/no**2, 1/no**2, 1/ne**2])
    # n = np.outer(n, n)

    n = np.array([1/2.40672491**2, 1/2.41373196**2, 1/2.38788325**2])
    n = np.outer(n, n)

    r = 0
    for m in range(len(freqs)):
        alpha_prime = alphas[m] * n
        #print(np.einsum('ij,k->ijk', alpha_prime, polarities[m]) / freqs[m]**2, '\n\n')
        r += np.einsum('ij,k->ijk', alpha_prime, polarities[m]) / freqs[m]**2
    
    r = r * 965060.241      # 1 e/(u A THz^2) = 965060.241 pm/V.
    r = -4 * np.pi * r
    
    r_voigt = np.zeros((6,3))
    r_voigt[0,:] = r[0,0,:]
    r_voigt[1,:] = r[1,1,:]
    r_voigt[2,:] = r[2,2,:]
    r_voigt[3,:] = r[1,2,:]
    r_voigt[4,:] = r[0,2,:]
    r_voigt[5,:] = r[0,1,:]

    return r, r_voigt



if __name__ == "__main__":

    alphas = []
    polarities = []
    freqs = []
    # vol = 12.4840002059999993 * 12.4840002059999993 * 3.9742000103000001

    born = get_Born(vasprun_born="../../03-Dielectric/KNL2/vasprun.xml")
    d_chi_d_tau = susceptibility_derivatives(0.01, 45)

    for band in np.arange(4, 10):
        
        freq, eig = get_phonons(qpoint_file="../../02-Phonons/cart/qpoints.yaml", band=band, poscar="POSCAR")
        polarity = calculate_polarity(born, eig)
        freqs.append(freq)
        polarities.append(polarity)

        alpha = calculate_Raman(d_chi_d_tau, eig)   
        alphas.append(alpha)

    freqs[0] = 1.19586
    r, r_voigt = calculate_Pockels_ionic(alphas, polarities, freqs)

    print(r_voigt)
    np.savetxt("IonicPockels.csv", r_voigt, delimiter=',')
    
    sns.heatmap(r_voigt, annot=True, center=0)
    plt.tick_params(labelsize=10, labelbottom = False, bottom=False, top=False, labeltop=True, left=False)
    plt.xticks(np.arange(3)+0.5, np.arange(1,4))
    plt.yticks(np.arange(6)+0.5, np.arange(1,7), rotation=0)
    plt.title(r"$r^{ion} (pm/V)$")
    plt.savefig("IonicPockels.pdf")
