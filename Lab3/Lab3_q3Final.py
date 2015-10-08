#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy

"""
Calculates the electric field from potential, and graphs the field Lines
over the equipotential surface 
"""

__author__ = "Eric Yeung"

Q = 1e-13 # coulombs
l = 1 # in millimeters
N = 10
from scipy.constants import epsilon_0 as e0

def V(r,z):
    
    def f(u):
        return exp(-(tan(u))**2)/((cos(u))**2*sqrt((z-l*tan(u))**2+r**2))
    
    a = -pi/2
    b = pi/2
    h = (b-a)/N

    s = f(a) + f(b) # from 5.9
    for i in range(1, N, 2):  # odd terms
        s += 4*f(a+i*h)
    for i in range(2, N-1, 2): # even terms
        s += 2*f(a+i*h)
    
    return Q/(4*pi*e0)*(h/3)*s

r = np.linspace(-6, 6, N)
z = np.linspace(-6, 6, N)

delta = 1e-15


gradr = np.empty([N, N])
gradz = np.empty([N, N])

gradr[0] = (V(r + delta, z) - V(r,z))/delta
gradr[-1] = (V(r,z) - V(r - delta, z))/delta
gradr[1:-2] = (V(r + delta, z) - V(r - delta,z))/(2*delta)

gradz[0] = (V(r, z + delta) - V(r,z))/delta
gradz[-1] = (V(r,z) - V(r, z - delta))/delta
gradz[1:-2] = (V(r, z + delta) - V(r,z - delta))/(2*delta)

fieldr = -1*gradr
fieldz = -1*gradz

R,Z = np.meshgrid(r, z)

plt.axes([-1.2, 1.2, -1.2, 1.2]) # Size of the plot
vecspace = 1

fill = plt.contourf(Z, R, V(R,Z), alpha=.8, cmap = 'Spectral')
C = plt.contour(Z, R, V(R,Z), colors='black')
plt.clabel(C, inline=True, fmt = '%.4f', fontsize=10)

q = plt.quiver(R[::vecspace], Z[::vecspace], fieldr[::vecspace], fieldz[::vecspace], color = 'deepskyblue', label = 'electric field')
plt.colorbar(fill,fraction=0.0338, pad=0.06)

plt.xlabel('z')
plt.ylabel('r')
plt.title('Electric Field Lines on Equipotential Plot')
plt.show()