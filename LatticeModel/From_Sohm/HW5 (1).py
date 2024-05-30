import numpy as np
import matplotlib.pyplot as plt


a = 5.43
b = 2*np.pi/a
# In units of 2*pi/a:
b1 = np.array([1,1,-1])*b
b2 = np.array([-1,1,1])*b
b3 = np.array([1,-1,1])*b

G0 = [0*b1]
G3 = [b1, -b1, b2, -b2, b3, -b3, b1+b2+b3, -b1-b2-b3]
G4 = [b1+b2, -b1-b2, b2+b3, -b2-b3, b3+b1, -b3-b1]
G8 = [b1-b2, b2-b1, b2-b3, b3-b2, b3-b1, b1-b3, 2*b1+b2+b3, -2*b1-b2-b3, b1+2*b2+b3, -b1-2*b2-b3, b1+b2+2*b3, -b1-b2-2*b3]

G = G0+G3+G4+G8
G = np.array(G)


# In units of 2*pi/a:
Gamma = np.array([0,0,0])*b
L = np.array([1/2,1/2,1/2])*b
X = np.array([0,1,0])*b

N = 50      # Number of points along each path.

# Path from L to Gamma:
delta = Gamma - L
path = np.outer(np.linspace(0,1,N,endpoint=False), delta) + L
# Path from Gamma to X:
delta = X - Gamma
path = np.vstack((path, np.outer(np.linspace(0,1,N), delta) + Gamma))
np.savetxt("k_new.txt", path.T)

E = [np.sum((G+k)**2, axis=1) for k in path]
E = np.array(E)


x = np.arange(2*N)
plt.plot(x,E)
plt.axvline(x = 0, label = 'L')
plt.axvline(x = N, label = 'Gamma')
plt.axvline(x = 2*N-1, label = 'X')
plt.xticks([0, N, 2*N-1], ['L', 'Gamma', 'X'])
plt.ylabel("Energy")
plt.title("Problem 1")
plt.savefig("HW5_prob1.pdf")
