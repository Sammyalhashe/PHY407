#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""Wave intensity after diffraction"""

__author__  = "Eric Yeung"

# Modified gaussxw.py for t
def gausstw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    t = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(t)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*t*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-t*p1)/(1-t*t)
        dt = p1/dp
        t -= dt
        delta = max(abs(dt))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-t*t)*dp*dp)

    return t, w

def C(t):
    return np.cos(0.5*np.pi*t*82)

def S(t):
    return np.sin(0.5*np.pi*t**2)

N = 50
a = 0

x = np.linspace(-5, 5)
b = np.sqrt(2/3)*x

#b = map(lambda x: np.sqrt(2/3)*x, x), lambda functions not so useful this time
#print b

t, w = gausstw(N)
tp = 0.5*(b-a)*t + 0.5*(b+a)
wp = 0.5*(b-a)*w

# Compute the first integral

s1 = []

for i in range(N):
    s1.append(wp[i]*C(tp[i]))

# Compute the second integral

s2 = []

for j in range(N):
    s2.append(wp[j]*S(tp[j]))


Iratio = [] 

for k in range(N):
    Iratio.append(1/8.*((2*s1[k] + 1)**2 + (2*s2[k] + 1)**2))

#print Iratio    

plt.plot(x, Iratio, color = 'r')
plt.xlabel('x')
plt.xlim([x.min(), x.max()])

plt.ylabel('$I/I_0$')
plt.title('Diffraction of Sound Wave')

plt.show()