#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt  
import numpy as np
from pylab import *

Q = 1e-13 # coulombs
l = 1 # in millimeters
N = 10
from scipy.constants import epsilon_0 as e0
from scipy import special as sp

def V(r,z):
    return Q/(4*pi*e0*l)*exp((r**2)/(2*l**2))*sp.kn(0, r**2/(2*l**2))

r = np.linspace(-0.25, 5, N)
z = np.linspace(-5, 5, N)

delta = 1e-5 # This isn't the same as vec_spacing
vec_spacing = 1

gradientr = np.zeros(N)
gradientz = np.zeros(N)

gradientr[0] = (V(r + delta, z) - V(r,z))/delta
gradientr[-1] = (V(r,z) - V(r - delta, z))/delta
gradientr[1:-2] = (V(r + delta, z) - V(r - delta,z))/(2*delta)

gradientz[0] = (V(r, z + delta) - V(r,z))/delta
gradientz[-1] = (V(r,z) - V(r, z - delta))/delta
gradientz[1:-2] = (V(r, z + delta) - V(r,z - delta))/(2*delta)

fieldr = -1*gradientr
fieldz = -1*gradientz


q = plt.quiver(r, z, fieldr, fieldz, color = 'dodgerblue', label = 'electric field')

plt.xlabel('r')
plt.ylabel('z')
plt.title('Electric Field')

plt.show()