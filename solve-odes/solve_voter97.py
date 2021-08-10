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
    dVdy = param[0]*np.cos(2*pi*x[0]) + 2*param[1]*pi*x[1]
     
    xDot[0] = x[2]
    xDot[1] = x[3]
    xDot[2] = -dVdx 
    xDot[3] = -dVdy
    
    return xDot


def py_momenta(x, y, px, param, TOTAL_ENERGY):
    
    potential_energy = pe_voter97(x, y, param)
    
    if TOTAL_ENERGY > (potential_energy + px**2/2):
        py = np.sqrt(2*(TOTAL_ENERGY - (potential_energy + px**2/2)))
    else:
        py = np.Nan
    
    return py

def px_momenta(x, y, py, param, TOTAL_ENERGY):
    
    potential_energy = pe_voter97(x, y, param)
    
    if TOTAL_ENERGY > (potential_energy + py**2/2):
        px = np.sqrt(2*(TOTAL_ENERGY - (potential_energy + py**2/2)))
    else:
        px = np.Nan
    
    return px


def pe_voter97(x, y, param):
    
    return np.cos(2*pi*x)*(1 + param[0]*y) + param[1]*pi*y**2
    

# Model parameters
d_1 = 4 #Coupling term
d_2 = 2*pi

param = np.array([d_1, d_2])

x_saddle = [1, -param[0]/(2*param[1]*pi)]
x_centre = [0.5, param[0]/(2*param[1]*pi)]

deltaE = 0.01

TOTAL_ENERGY = pe_voter97(x_saddle[0],x_saddle[1], param) + deltaE

y_shift = 0.1

n=5
N=n**2


# Define initial conditions and solver parameters
C=84
omega = np.sqrt(2*param[1]*pi)
T= C*(2*pi/omega)
t_span = [0, T]

x0_vals = np.repeat(x_centre[0], n)
y0_vals = np.linspace(x_centre[1]-y_shift, x_centre[1] + y_shift, num=n)
py0_vals = np.linspace(-1, 1, num=n)
y0_mesh,py0_mesh = np.meshgrid(y0_vals,py0_vals)

px0 = lambda y,py: np.sqrt(2*(TOTAL_ENERGY - (pe_voter97(x_centre[0], y, param) + (py**2)/2)))            
px0_vals =  px0(y0_mesh, py0_mesh)

data = np.empy(n,4)


# positive momenta to cross the saddle
 



#init_conds = [x0, y0, px0, py0]


relTol = 1e-12
absTol = 1e-9

#sol = solve_ivp(ham2dof_voter97, t_span, init_conds, method='RK45', rtol = relTol, atol = absTol, args = param)


fig, ax = plt.subplots()

def traj2d(x, y, x0, y0, param):
    plt.plot(x, y, '-b')
    plt.plot(x0,y0,'gx')
    plt.plot(x[-1],y[-1],'rx')
    X = np.arange(0, 3, 0.001)
    Y = np.arange(-1, 1, 0.001)
    X, Y = np.meshgrid(X, Y)
    Z = np.cos(2*pi*X)*(1+d_1*Y) + d_2*pi*Y**2
    contour = plt.contourf(X, Y, Z, 20, cmap='hot')
    ax.set_xlabel(r'$q_{1}$')
    ax.set_ylabel(r'$q_{2}$')
    plt.colorbar()
    plt.show()


def traj3d(x,y,x0,y0,param):
    ax = plt.axes(projection='3d')
    ax.set_xlabel(r'$q_{1}$')
    ax.set_ylabel(r'$q_{2}$')
    ax.set_zlabel(r'$V(q_{1},q_{2})$')
    ax.plot3D(x, y, pe_voter97(x,y,param), 'black')
    ax.scatter(x0,y0,pe_voter97(x0,y0,param), c='green')
    ax.scatter(x[-1],y[-1], pe_voter97(x[-1], y[-1], param), c='red')
    
#X = np.arange(min(sol.y[0,:]), max(sol.y[0,:]), 0.001)
#Y = np.arange(min(sol.y[1,:]), max(sol.y[1,:]), 0.001)
x,y=sol.y[0,:],sol.y[1,:]

#traj2d(x, y, init_conds[0], init_conds[1], param)

#traj3d(x, y, init_conds[0], init_conds[1], param)

def intermediate_region_crossed(x, y, param, R):
    if ((x-3/2)**2 + (y-param[0]/(2*param[1]*pi))**2).any() <= R**2:
        print('yeh')
    else:
        print('nah')
        
def intermediate_region_crossed2(x, y, param, R):
    for i in range(n):
        if x[i] in ():
            print('yeh')
    else:
        print('nah')        

#ax.plot_wireframe(X,Y,Z,rstride=5, cstride=5)
intermediate_region_crossed(x,y, param, 0.1)


#ax.plot3D(sol.y[0,:], sol.y[1,:],pe_voter97(sol.y[0,:],sol.y[1,:],param), '-b')

#ax = plt.gca()
#ax.axes.xaxis.set_ticklabels([])
#ax.axes.yaxis.set_ticklabels([])
#ax.axes.zaxis.set_ticklabels([])

#ax.xaxis.labelpad=-10
#ax.yaxis.labelpad=-10
#ax.zaxis.labelpad=-10