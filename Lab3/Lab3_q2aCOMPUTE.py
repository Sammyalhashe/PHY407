#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from gaussxw import gaussxw

"""
Compares period to classical limit 
"""

__author__ = "Eric Yeung"

def g(x):
    return 1/(c*np.sqrt(abs(0.5*k*(x0**2 - x**2)*(2*m*c**2 + 0.5*k*(x0**2 - x**2))/(m*c**2 + 0.5*k*(x0**2 - x**2)))))

# From the first assignment
k = 12.
m = 1. 
c = 2.998e8
xc = 8.654e7

N = 8
a = 0.0
x0 = 0.01
#x0 = np.linspace(0.01, 1, N)
#x0 = np.linspace(1, 10*xc, N)
b = x0

x,w = gaussxw(N)
xp = 0.5*(b-a)*x + 0.5*(b+a)
wp = 0.5*(b-a)*w

T = 0
for k in range(1,N):
    T += wp[k]*g(xp[k])

"""
T = wp*g(xp)
"""
    
period = 4*T; print period

classicalLimit = 2*np.pi*np.sqrt(m/k); print classicalLimit

relativisticLimit = 4*x0/c; print relativisticLimit


plt.plot(w*g(x), g(x))
plt.xlabel('wg(x)')
plt.ylabel('g(x)')

plt.show()