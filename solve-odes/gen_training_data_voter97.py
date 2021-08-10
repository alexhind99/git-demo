# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 12:30:21 2021

@author: alexa
"""
import voter97_ham as voter
from scipy.integrate import solve_ivp
import numpy as np
from math import pi as pi 
import matplotlib.pyplot as plt


# Model parameters(Optimising Transition States paper)
d_1 = 4 #Coupling term
d_2 = 2*pi

param = np.array([d_1, d_2])

# Equilibrium points
x_saddle = [1, -param[0]/(2*param[1]*pi)]
x_centre = [0.5, param[0]/(2*param[1]*pi)]


deltaE = 0.1
TOTAL_ENERGY = voter.pe_voter97(x_saddle[0],x_saddle[1], param) + deltaE  #Define initial conditions and solver parameters

#Integrate for T
#stop if trajectory reaches region of intermediate saddle. 
C=10
omega = np.sqrt(2*param[1]*pi)
T= C*(2*pi/omega)
t_span = [0, T]



#Generate training data
n = 100
N=n**2  #Sample size
data = np.zeros((N,4))

y_shift = 0.2
delta = 0.1

x0_vals = np.repeat(x_centre[0], n)
y0_vals = np.linspace(x_centre[1] - y_shift, x_centre[1] + 0.1, num=n)
py0_vals = np.linspace(-1.5, 0.5, num=n)
y0_mesh, py0_mesh = np.meshgrid(y0_vals,py0_vals)

px0 = lambda y,py: np.sqrt(2*(TOTAL_ENERGY - (voter.pe_voter97(x_centre[0], y, param) + (py**2)/2)))            
px0_vals =  px0(y0_mesh, py0_mesh)




data[:,0] = np.repeat(x0_vals,n)
data[:,1] = np.tile(y0_vals,n)
data[:,2] = px0_vals.flatten()
data[:,3] = np.repeat(py0_vals,n)

relTol = 1e-12
absTol = 1e-9

def event_escape_left(t, states, *params):
    return states[0] - (-1.25)
event_escape_left.terminal = True
event_escape_left.direction = 0

def region_reached(t, x, param):
    return x[0]-1

region_reached.terminal = True
region_reached.direction = 0

def intermediate_well_event(t, x, *param):
    return (voter.intermediate_region_crossed(x[0],x[1],param, 0.5))

intermediate_well_event.terminal = True
intermediate_well_event.terminal = 0


sol = solve_ivp(voter.ham2dof_voter97, t_span, data[150], method='RK45', \
                rtol = relTol, atol = absTol, args = param, dense_output = True, event = lambda t, x: region_reached(t, x, *param))
                     
x, y = sol.y[0,:], sol.y[1,:]
                        
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


traj2d(sol.y[0,:], sol.y[1,:], data[150,0], data[150,1], param)



def gen_train_data(data, param, t_span):
    react_data = np.empty(len(data))
    for i in range(len(data)):
        sol = solve_ivp(voter.ham2dof_voter97, t_span, data[i], method='RK45', \
                rtol = relTol, atol = absTol, args = param)
        #region = np.where(voter.intermediate_region_crossed(sol.y[0,:],sol.y[1,:],param, 0.1), \
        #                  voter.intermediate_region_crossed(sol.y[0,:],sol.y[1,:],param, 0.1) > 0, True)
        saddle0 = np.where(sol.y[0,:], sol.y[0,:] < 0, True)
        saddle1 = np.where(sol.y[0,:], sol.y[0,:] > 1, True)  
        if np.any(saddle0):
            react_data[i] = 2
        if np.any(saddle1):
            react_data[i] = 1
        else:
            react_data[i] = 0
    
    #np.concatenate((react_data.T, data), axis=1)
    return react_data


#react_data = gen_train_data(data, param, t_span)

#ax.scatter(data[:,1], data[:,3], c=react_data, cmap='inferno', s=20, edgecolors='k')
#ax.set_xlabel(r'$q_{2}$')
#ax.set_ylabel(r'$p_{2}$')
plt.show()



    



     