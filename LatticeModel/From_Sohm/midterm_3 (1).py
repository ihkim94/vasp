import numpy as np
from numpy import cos, sin, pi
import matplotlib.pyplot as plt


# Tight binding parameters in eV.
Es=0; Ep=7.2
Vss=-8.13; Vsp=5.88; Vxx=3.17; Vxy=7.51     # Vxx=1.71?


def E(k):
    # In units of 2*pi/a.
    g0 = +cos(pi*k[0]/2) * cos(pi*k[1]/2) * cos(pi*k[2]/2) -1j * sin(pi*k[0]/2) * sin(pi*k[1]/2) * sin(pi*k[2]/2)
    g1 = -cos(pi*k[0]/2) * sin(pi*k[1]/2) * sin(pi*k[2]/2) +1j * sin(pi*k[0]/2) * cos(pi*k[1]/2) * cos(pi*k[2]/2)
    g2 = -sin(pi*k[0]/2) * cos(pi*k[1]/2) * sin(pi*k[2]/2) +1j * cos(pi*k[0]/2) * sin(pi*k[1]/2) * cos(pi*k[2]/2)
    g3 = -sin(pi*k[0]/2) * sin(pi*k[1]/2) * cos(pi*k[2]/2) +1j * cos(pi*k[0]/2) * cos(pi*k[1]/2) * sin(pi*k[2]/2)

    H = np.zeros((8,8), dtype=complex)
    H[0,0] = H[1,1] = Es
    H[2,2] = H[3,3] = H[4,4] = H[5,5] = H[6,6] = H[7,7] = Ep
    H[0,1] = Vss*g0
    H[0,5] = Vsp*g1
    H[0,6] = Vsp*g2
    H[0,7] = Vsp*g3
    H[1,2] = -Vsp*np.conjugate(g1)
    H[1,3] = -Vsp*np.conjugate(g2)
    H[1,4] = -Vsp*np.conjugate(g3)
    H[2,5] = Vxx*g0
    H[2,6] = Vxy*g3
    H[2,7] = Vxy*g1
    H[3,5] = Vxy*g3
    H[3,6] = Vxx*g0
    H[3,7] = Vxy*g1
    H[4,5] = Vxy*g1
    H[4,6] = Vxy*g2
    H[4,7] = Vxx*g0

    for i in range(8):
        for j in range(i+1,8):
            H[j,i] = np.conjugate(H[i,j])

    return np.linalg.eigvals(H)



# In units of 2*pi/a:
Gamma = np.array([0,0,0])
L = np.array([1/2,1/2,1/2])
X = np.array([1,0,0])

N = 100      # Number of points along each path.

# Path from L to Gamma:
delta = Gamma - L
path = np.outer(np.linspace(0,1,N,endpoint=False), delta) + L
# Path from Gamma to X:
delta = X - Gamma
path = np.vstack((path, np.outer(np.linspace(0,1,N), delta) + Gamma))

energies = [E(k) for k in path]
energies = np.vstack(energies)
energies = np.sort(energies, axis=1)


# Effective mass:
cond = energies[:,4]
min = np.argmin(cond)   # conduction band minimum
k_min = path[min]

delta = X - Gamma
h = delta/1000
E_min = np.real(np.sort(E(k_min))[4])
E_min_p = np.real(np.sort(E(k_min+h))[4])
E_min_m = np.real(np.sort(E(k_min-h))[4])

d2E = (E_min_p - 2*E_min + E_min_m) / np.linalg.norm(h)**2     # finite difference second derivative
m = (0.1973*11571)**2 / d2E   # The coefficient converts the units of m to eV.
m = m/0.51E+6   # Dividing by the rest mass of an electron
print(m)


# Band Structure plot:
x = np.arange(2*N)
plt.plot(x, energies)
plt.plot(x[min],E_min, 'o')
plt.axvline(x = 0, label = 'L')
plt.axvline(x = N, label = 'Gamma')
plt.axvline(x = 2*N-1, label = 'X')
plt.xticks([0, N, 2*N-1], ['L', 'Gamma', 'X'])
plt.ylabel("Energy")
plt.show()
