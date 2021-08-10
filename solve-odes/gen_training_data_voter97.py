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
import time

start_time = time.time()

# Model parameters(Optimising Transition States paper 4, 2*pi)
d_1 = 0 #Coupling term
d_2 = 1

param = np.array([d_1, d_2])

# Equilibrium points
x_saddle = [1, -param[0]/(2*param[1]*pi)]
x_centre = [0.5, param[0]/(2*param[1]*pi)]


deltaE = 0.1
TOTAL_ENERGY = voter.pe_voter97(x_saddle[0],x_saddle[1], param) + deltaE  #Define initial conditions and solver parameters

#Integrate for T
#stop if trajectory reaches region of intermediate saddle. 
C=5
omega = np.sqrt(2*param[1]*pi)
T= C*(2*pi/omega)
t_span = [0, T]



#Generate training data
n = 50
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



def saddle1_crossed(t, x, *param):
    return x[0]-1

saddle1_crossed.terminal = True
saddle1_crossed.direction = 0

def intermediate_well_event(t, x, *param):
    return (voter.intermediate_region_crossed(x[0],x[1],param, 0.1))

intermediate_well_event.terminal = True
intermediate_well_event.terminal = -1    #Exterior to interior


#sol = solve_ivp(voter.ham2dof_voter97, t_span, data[150], method='RK45', \
 #               rtol = relTol, atol = absTol, args = param, dense_output = True, events = (intermediate_well_event))
                     
#x, y = sol.y[0,:], sol.y[1,:]
                        
fig, ax = plt.subplots()


#voter.traj2d(sol.y[0,:], sol.y[1,:], data[150,0], data[150,1], param)



def gen_train_data(data, param, t_span):
    react_data = np.empty(len(data))
    for i in range(len(data)):
        sol = solve_ivp(voter.ham2dof_voter97, t_span, data[i], method='RK45', \
                rtol = relTol, atol = absTol, args = param,  dense_output = True, events = (intermediate_well_event))
        saddle0 = np.where(sol.y[0,:], sol.y[0,:] < 0, True)
        
        if np.any(saddle0):
            react_data[i] = np.NaN
        if sol.status == 1:
            react_data[i] = 1
        else:
            react_data[i] = 0

    #np.concatenate((react_data.T, data), axis=1)
    return react_data


react_data = gen_train_data(data, param, t_span)

ax.scatter(data[:,1], data[:,3], c=react_data, cmap='inferno', s=20, edgecolors='k')
ax.set_xlabel(r'$q_{2}$')
ax.set_ylabel(r'$p_{2}$')
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

    



     