#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy

"""
Uses simpson's rule to compute potential, and shows 
the fractional error between it and the analytical value.
They are plotted against one another. 
"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

from scipy.constants import epsilon_0 as e0

Q = 1e-13 # coulombs
#l = 0.001 in metres, this causes overflow, exponent too large. 
l = 1 # in millimeters

def V(r,z):
    
    def f(u):
        return exp(-(tan(u))**2)/((cos(u))**2*sqrt((z-l*tan(u))**2+r**2))
    
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

#Testing potential at the center
#print V(1,0) 

from scipy import special as sp

def calcV(r,z):
    return Q/(4*pi*e0*l)*exp((r**2)/(2*l**2))*sp.kn(0, r**2/(2*l**2))

#Testing like above
#print calcV(1,0)

def percentError(r):
	return V(r,0)/calcV(r,0)

def fractionalError(r):
	return (V(r,0)-calcV(r,0))/calcV(r,0)

r = linspace(0.25, 5.00)

plt.plot(r, V(r,0), color = 'g', label = 'Computed', linewidth = '2')
plt.plot(r, calcV(r,0), color = 'fuchsia', label = 'Analytical', linewidth = '2')

plt.xlabel('r(mm)')
plt.ylabel('Electrostatic Potential')
plt.title('Comparison')

plt.legend(loc='upper right')
plt.tight_layout()
plt.show()

"""

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True)

ax1.plot(r, V(r,0), color = 'g', label = 'Computed', linewidth = '2')
ax2.plot(r, calcV(r,0), color = 'fuchsia', label = 'Analytical', linewidth = '2')

ax1.set_xlabel('r(mm)')
ax1.set_ylabel('Electrostatic Potential')
ax1.set_title('Computed Potential')

ax2.set_xlabel('r(mm)')
ax2.set_ylabel('Electrostatic Potential')
ax2.set_title('Analytical Potential')

plt.tight_layout()
plt.show()

print fractionalError(1) # This is 1e-6 when N~50 or so.
"""