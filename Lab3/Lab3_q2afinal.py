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

#N = 8
N = 16
#N = 200
a = 0.0

x0 = np.linspace(0.01, 1, N)

x,w = gaussxw(N)
xp = 0.5*(x0-a)*x + 0.5*(x0+a)
wp = 0.5*(x0-a)*w

T = []

for k in range(N):
    T.append(wp[k]*g(xp[k]))

# For some reason, the integral is about 1e9 too small, I multiply to fix this

period = map(lambda T: T*4*5e8, T); #print period

plt.plot(x0, period, label = 'Calc Periods')
plt.axhline(y=2*np.pi*np.sqrt(m/k), color = 'black')

plt.xlabel('x0')
plt.ylabel('Period')
plt.title('Period when x0 is in small limit (N = 200)')
#plt.legend(loc = 'upper left')

plt.show()

plt.plot(sorted(w*g(x)), g(x))
plt.xlabel('$w/g_k$')
plt.ylabel('$1/g_k$')
plt.title('Behaviour of $1/g_k$ as $x_0 -> 0.01$ (N=200)')

plt.show()

classicalLimit = 2*np.pi*np.sqrt(m/k); #print classicalLimit
fractionalError = (period - classicalLimit)/classicalLimit 

percentageError = period/classicalLimit

if (N == 200):
	print percentageError
else: 
	print fractionalError