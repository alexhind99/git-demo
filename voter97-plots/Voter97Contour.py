import matplotlib.pyplot as plt
import numpy as np
from math import pi as pi

d1=4
d2=2*pi

X = np.arange(0, 3, 0.001)
Y = np.arange(-1, 1, 0.001)
X, Y = np.meshgrid(X, Y)
Z = np.cos(2*pi*X)*(1+d1*Y) + d2*pi*Y**2

fig, ax = plt.subplots()
contour = plt.contourf(X, Y, Z, 50, cmap='hot')

ax.set_xlabel(r'$q_{1}$')
ax.set_ylabel(r'$q_{2}$')

plt.colorbar()
plt.show()
