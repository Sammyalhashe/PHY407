#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

u = np.linspace(0.980, 0.984, 10000) # 10000 even spaces

def p(u):
    return (1-u)**8

def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

plt.plot(u, abs(p(u)-q(u))/abs(p(u)), color = 'dodgerblue')
plt.xlabel('u')
plt.xlim([u.min(), u.max()])

plt.ylabel('abs(p-q)/abs(p)')
plt.title('Divergence of Fractional Error as u $\longrightarrow 1$')

plt.show()

"""
Old Print method (Ignore)

u = 0.980

while u <= 0.984:
    p = (1-u)**8
    q = 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8
    print u, abs(p-q)/abs(p)
    u += 0.0001

"""