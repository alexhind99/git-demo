# -*- coding: utf-8 -*-
"""
Solving Hamilton's equations numerically for the Voter97 Potential
"""
import numpy as np
from scipy.integrate import solve_ivp
from math import pi as pi
import matplotlib.pyplot as plt   


def ham2dof_voter97(t, x, *param):
    
    xDot = np.zeros((4,))
    
    dVdx = -2*pi*np.sin(2*pi*x[0])*(1 + param[0]*x[1])
    dVdy = param[0]*np.cos(2*pi*x[0]) + 2*param[1]*x[1]
        
    xDot[0] = x[2]
    xDot[1] = x[3]
    xDot[2] = -dVdx 
    xDot[3] = -dVdy
    
    return xDot


def py_momenta(x, y, px, param):
    
    potential_energy = pe_voter97(x, y, param)
    
    if param[2] > (potential_energy + px**2/2):
        py = np.sqrt(2*(param[2]- (potential_energy + px**2/2)))
    else:
        py = np.Nan
    
    return py


def pe_voter97(x, y, param):
    
    return np.cos(2*pi*x)*(1 + param[0]*y) + param[1]*pi*y**2
    
#and (py**2 > 2*total - 4*d_2*pi*y**2)

def reactive_uncoup(x, y, px, py, total, d_2):
    if (abs(px) > 2*np.sin(pi*x)):
        print('R')
    else:
        print('NR')

# Model parameters
d_1 = 4
d_2 = 1

TOTAL_ENERGY = 3
param = np.array([d_1, d_2, TOTAL_ENERGY])


# Define initial conditions and solver parameters
t_span = [0, 10]
x0 = 0.5
y0 = -0.2
px0 = -1.9
py0 = py_momenta(x0, y0, px0, param)# positive momenta to cross the saddle

init_conds = [x0, y0, px0, py0]


relTol = 1e-12
absTol = 1e-9

sol = solve_ivp(ham2dof_voter97, t_span, init_conds, method='RK45', rtol = relTol, atol = absTol, args = param)



reactive_uncoup(x0, y0, px0, py0, TOTAL_ENERGY, d_2)

fig, ax = plt.subplots()
plt.plot(sol.y[0,:], sol.y[1,:], '-b')

plt.plot(x0,y0,'gx')
plt.plot(sol.y[0,-1],sol.y[1,-1],'mx')

X = np.arange(min(sol.y[0,:]), max(sol.y[0,:]), 0.001)
Y = np.arange(min(sol.y[1,:]), max(sol.y[1,:]), 0.001)
X, Y = np.meshgrid(X, Y)
Z = np.cos(2*pi*X)*(1+d_1*Y) + d_2*pi*Y**2

contour = plt.contourf(X, Y, Z, 20, cmap='hot')

ax.set_xlabel(r'$q_{1}$')
ax.set_ylabel(r'$q_{2}$')

plt.colorbar()
plt.show()

#fig.savefig('reactivetraj2DS.png', dpi=900)
