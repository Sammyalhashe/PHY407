#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt  
import numpy as np
from pylab import *
import scipy

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

r = np.linspace(-0.25, 5, N)
z = np.linspace(-5, 5, N)

delta = 1e-5

R,Z = np.meshgrid(r, z)

plt.axes([-1, 1, -1, 1]) # Size of the plot

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

q = plt.quiver(R, Z, fieldr, fieldz, color = 'dodgerblue', label = 'electric field')
# vector spacing different from delta

plt.xlabel('r')
plt.ylabel('z')
plt.title('Electric Field')

plt.show()


plt.quiver = (R,Z, V(r,z))

plt.show()
