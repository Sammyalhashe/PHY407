#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy

"""
Script to contour plot the electrostatic potential
"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

Q = 1e-13 # coulombs
l = 1 # in millimeters

from scipy.constants import epsilon_0 as e0

# From part A

def V(r,z):
    
    def f(u):
        return exp(-(tan(u))**2)/((cos(u))**2*np.sqrt((z-l*tan(u))**2+r**2))
    
    N = 8
    a = -pi/2
    b = pi/2
    h = (b-a)/N

    s = f(a) + f(b) # from 5.9
    for i in range(1, N, 2):  # odd terms
        s += 4*f(a+i*h)
    for i in range(2, N-1, 2): # even terms
        s += 2*f(a+i*h)
    
    return Q/(4*pi*e0)*(h/3)*s

N= 8

r = np.linspace(-0.25, 5, N)
z = np.linspace(-5, 5, N)
R,Z = np.meshgrid(r, z)

plt.axes([-1.5, 1.5, -1.5, 1.5]) # Size of the plot

plt.contourf(Z, R, V(R, Z), alpha=.75, cmap = 'Spectral')
C = plt.contour(Z, R, V(R, Z), colors='black')
plt.clabel(C, inline=True, fmt = '%.4f', fontsize=10)

plt.xlabel('z')
plt.tick_params(axis='y', labelleft='off', labelright = 'on')
plt.ylabel('r')
plt.title('Potential of Line Charge')

plt.show()