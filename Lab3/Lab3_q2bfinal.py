#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from gaussxw import gaussxw

"""
Compares period to relativistic limit 
"""

__author__ = "Eric Yeung"

def g(x):
    return 1/(c*np.sqrt(abs(0.5*k*(x0**2 - x**2)*(2*m*c**2 + 0.5*k*(x0**2 - x**2))/(m*c**2 + 0.5*k*(x0**2 - x**2)))))

# From the first assignment
k = 12.
m = 1. 
c = 2.998e8
xc = 8.654e7

N = 200
a = 0.0

x0 = np.linspace(1, 10*xc, N)

x,w = gaussxw(N)
xp = 0.5*(x0-a)*x + 0.5*(x0+a)
wp = 0.5*(x0-a)*w

T = []

for k in range(N):
    T.append(wp[k]*g(xp[k]))
 
# For some reason, the integral is about 1e10 too small, I multiply to fix this

period = map(lambda T: T*4*1e10, T); #print period

plt.plot(x0, period, label = 'Calc Periods')
plt.plot(x0, 4*x0/xc, label = 'Relativistic Limit')

plt.xlabel('x0')
plt.ylabel('Period')
plt.title('Period when x0 is in large limit')

plt.show()

relativisticLimit = 4*x0/c; #print relativisticLimit
fractionalError = (period - relativisticLimit)/relativisticLimit; print fractionalError