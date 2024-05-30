import numpy as np
import matplotlib.pyplot as plt

# Problem 1:
a=0; b=-1

E = []
for n in range(3,11):
    H = np.zeros((n,n))
    for i in range(n):
        H[i,i] = a
        H[i,(i+1)%n] = b
        H[i,(i-1)%n] = b

    E.append(np.linalg.eigvals(H))


print(E)

x=[]; y=[]
for i in range(len(E)):
    for j in range(len(E[i])):
        x.append(i+3)
        y.append(E[i][j])
plt.plot(x,y,'.')
plt.show()
