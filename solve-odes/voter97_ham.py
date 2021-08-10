# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 13:13:07 2021

@author: alexa
"""

import numpy as np
from scipy.integrate import solve_ivp
from math import pi as pi
import matplotlib.pyplot as plt  
 
def pe_voter97(x, y, param):
    
    return np.cos(2*pi*x)*(1 + param[0]*y) + param[1]*pi*y**2

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


def intermediate_region_crossed(x,y, param, R):
        return ((x-3/2)**2 + (y-param[0]/(2*param[1]*pi))**2) - R**2
        

def saddle0_cross(x):
    if x < 0:
        return True
    
def saddle1_cross(t,x):
    if x > 1:
        return 0
    else:
        return 1
    
def region_reached1(t,x, param):
    return x[0]-1
