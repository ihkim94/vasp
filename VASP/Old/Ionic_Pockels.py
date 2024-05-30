from pymatgen.io.vasp import Vasprun
import numpy as np

number_of_sites = 5
delta = 0.005 # From Generate_disp_Ionic_Pockels
dchi_dtau = np.zeros((3,3))

#####
# summation of dchi_dtau 
for i in range(0, number_of_sites):
    vrun_plus = Vasprun(f'1.x/1.plus/{i+1:02d}/vasprun.xml')
    vrun_minus = Vasprun(f'1.x/2.minus/{i+1:02d}/vasprun.xml')
    die_plus = vrun_plus.epsilon_static
    die_minus = vrun_minus.epsilon_static
    print(f'{i+1:02d} site:', (np.array(die_plus) - np.array(die_minus))/(2*delta))
    dchi_dtau += (np.array(die_plus) - np.array(die_minus))/(2*delta)
    print(dchi_dtau)

for i in range(0, number_of_sites):
    vrun_plus = Vasprun(f'2.y/1.plus/{i+1:02d}/vasprun.xml')
    vrun_minus = Vasprun(f'2.y/2.minus/{i+1:02d}/vasprun.xml')
    die_plus = vrun_plus.epsilon_static
    die_minus = vrun_minus.epsilon_static
    dchi_dtau += (np.array(die_plus) - np.array(die_minus))/(2*delta)

for i in range(0, number_of_sites):
    vrun_plus = Vasprun(f'3.z/1.plus/{i+1:02d}/vasprun.xml')
    vrun_minus = Vasprun(f'3.z/2.minus/{i+1:02d}/vasprun.xml')
    die_plus = vrun_plus.epsilon_static
    die_minus = vrun_minus.epsilon_static
    dchi_dtau += (np.array(die_plus) - np.array(die_minus))/(2*delta)
#####
print("")
print("")
print("Summation of dchi_dtau:")
print(dchi_dtau)
# print("Static dielectric tensor:")
# print((dielectric_tensor))
##################################
from pymatgen.io.vasp import Vasprun
import numpy as np

BORN_1 = np.array([[2.71521, 0.00000000, 0.00000000], [0.00000000, 2.71521, 0.00000000], [0.00000000, 0.00000000, 2.76603]])
BORN_2 = np.array([[7.56336, 0.00000000, 0.00000000], [0.00000000, 7.56336, 0.00000000], [0.00000000, 0.00000000, 6.63160]])
BORN_3 = np.array([[-6.05159, 0.00000000, 0.00000000], [0.00000000, -2.14690, 0.00000000], [0.00000000, 0.00000000, -2.01708]])
BORN_4 = np.array([[-2.14690, 0.00000000, 0.00000000], [0.00000000, -6.05159, 0.00000000], [0.00000000, 0.00000000, -2.01708]])  
BORN_5 = np.array([[-2.08007, 0.00000000, 0.00000000], [0.00000000, -2.08007, 0.00000000], [0.00000000, 0.00000000, -5.36346]])


eigv_1 = np.array([[ -0.76732911023526],
[  0.00000000009508],
[  0.00000000000991]])

eigv_2 = np.array([[ -0.45518878750146],
[  0.00000000005614],
[  0.00000000000586]])

eigv_3 = np.array([[ -0.26045086937023],
[  0.00000000003246],
[  0.00000000000339]])

eigv_4 = np.array([[ -0.26084599859961],
[  0.00000000003246],
[  0.00000000000339]])

eigv_5 = np.array([[ -0.26102473824789],
[  0.00000000003246],
[  0.00000000000339]])

polar_sum = np.dot(BORN_1, eigv_1) + np.dot(BORN_2,eigv_2) + np.dot(BORN_3,eigv_3) + np.dot(BORN_4,eigv_4) + np.dot(BORN_5,eigv_5)

print(polar_sum)
####################################
Thz_to_cm = 33.356
phonon_frequency = np.array[0.0082379220, 4.8253696242, 5.1979869319, 5.2329837917, 6.5059899262, 8.5066523784, 8.5087657160, 8.5243068262, 13.2640517300, 14.1353586741, 14.6818988090, 20.9128593710] # excluding first 3 imaginary modes
Ionic_pockels = 4 * np.pi * 1/(0.0082379220*Thz_to_cm)**2 * dchi_dtau * polar_sum
print(Ionic_pockels)