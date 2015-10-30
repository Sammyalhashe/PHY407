#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from gaussxw import gaussxw

def g(x):
    return 1/(c*np.sqrt(abs(0.5*k*(x0**2 - x**2)*(2*m*c**2 + 0.5*k*(x0**2 - x**2))/(m*c**2 + 0.5*k*(x0**2 - x**2)))))

# From the first assignment
k = 12.
m = 1. 
c = 2.998e8
xc = 8.654e7

#N = 200
N = 8
a = 0.0
# Let x0 be 1 cm
#x0 = np.linspace(1, 10*xc, N)
x0 = np.linspace(0.01, 1, N)

x,w = gaussxw(N)
xp = 0.5*(x0-a)*x + 0.5*(x0+a)
wp = 0.5*(x0-a)*w

T = []

for k in range(N):
    T.append(wp[k]*g(xp[k]))

#period = 4*T; print period
period = map(lambda T: T*4, T); print period

plt.plot(x0, period)

plt.xlabel('x0')
plt.ylabel('Period')

plt.title('Period when x0 is in small limit')
plt.show()

#plt.plot(x0,T,x0,2*np.pi*(m/k)**0.5+0*T, x0, 4*x0/c)

"""

T = 0.0
for k in range(N):
    T += wp[k]*g(xp[k])

print(T*4)

classicalLimit = 2*np.pi*np.sqrt(m/k); print classicalLimit

relativisticLimit = 4*x0/c; print relativisticLimit

"""
