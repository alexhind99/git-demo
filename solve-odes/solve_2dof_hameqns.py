#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical solution of odes (Hamilton's equations for 2 dof systems)
                            
"""

# Import scipy, matplotlib, numpy
import numpy as np
from scipy.integrate import solve_ivp


# Define vector field for the system
def ham2dof_deleonberne(t, x, *parameters):
    """ 
    Returns the Hamiltonian vector field (Hamilton's equations of motion) 
    
    Used for passing to ode solvers for integrating initial conditions over a time interval.
    Parameters
    ----------
    t : float
        time instant
    x : float (list of size 4)
        phase space coordinates at time instant t
    parameters : float (list)
        model parameters
    Returns
    -------
    xDot : float (list of size 4)
        right hand side of the vector field evaluated at the phase space coordinates, x, at time instant, t
    """
    
    xDot = np.zeros((4,))
    
    dVdx = -2*parameters[3]*parameters[4]*np.exp(-parameters[4]*x[0])*(np.exp(-parameters[4]*x[0]) - 1) - \
        4*parameters[5]*parameters[4]*x[1]**2*(x[1]**2 - 1)*np.exp(-parameters[5]*parameters[4]*x[0])
    dVdy = 8*x[1]*(2*x[1]**2 - 1)*np.exp(-parameters[5]*parameters[4]*x[0])
        
    xDot[0] = x[2]/parameters[0]
    xDot[1] = x[3]/parameters[1]
    xDot[2] = -dVdx 
    xDot[3] = -dVdy
    
    return xDot


def py_momenta(x, y, px, parameters):
    """
    Returns the y-momentum for a fixed energy given the other 3 coordinates
    """
    
    total_energy = parameters[-1]
    potential_energy = pe_deleonberne(x, y, parameters)
    if total_energy > (potential_energy + (px**2)/(2*parameters[0])):
        py = np.sqrt(2*parameters[1]*(total_energy \
                                  - (potential_energy + (px**2)/(2*parameters[0]))))
    else:
        py = np.Nan
            
    
    return py


def pe_deleonberne(x, y, parameters):
    """ Returns the potential energy at the configuration space coordinates 
    
    Parameters
    ----------
    x : float
        configuration space coordinate
    y : float
        configuration space coordinate
    parameters : float (list)
        model parameters
    Returns
    -------
    float 
        potential energy of the configuration
    
    """
    
    return parameters[3]*( 1 - np.exp(-parameters[4]*x) )**2 \
        + 4*y**2*(y**2 - 1)*np.exp(-parameters[5]*parameters[4]*x) \
        + parameters[2]

# Set parameters for the model system
MASS_A = 8.0 
MASS_B = 8.0 # De Leon, Marston (1989)
EPSILON_S = 1.0
D_X = 10.0
ALPHA = 1.00
LAMBDA = 1.5
TOTAL_ENERGY = 1.1
parameters = np.array([MASS_A, MASS_B, EPSILON_S, D_X, LAMBDA, ALPHA, TOTAL_ENERGY])

# Define initial conditions and solver parameters
t_span = [0, 50]
x0 = 0.1
y0 = -0.2
px0 = -0.1
py0 = py_momenta(x0, y0, px0, parameters)# positive momenta to cross the saddle
y0 = [x0, y0, px0, py0]

relTol = 1e-12
absTol = 1e-9

# Solve for trajectories
sol = solve_ivp(ham2dof_deleonberne, t_span, y0, method='RK45', \
                rtol = relTol, atol = absTol, args = parameters)

# Visualise on the 2D potential
import matplotlib.pyplot as plt                                                                                  

fig = plt.figure(figsize = (5,5))
plt.plot(sol.y[0,:], sol.y[1,:], '-r') 
# plt.scatter(sol.y[0,0], sol.y[1,0],'xg')
# plt.scatter(sol.y[0,-1], sol.y[1,-1], 'xr')

plt.show()

#%% 

 


