#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

"""Compares the noise in p(u) and q(u)"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

def p(u):
    return (1-u)**8


def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

N = 1000 # number of points
u = np.random.uniform(0.98,1.02,N)

f,(ax1, ax2) = plt.subplots(1, 2, sharex='row')
ax1.plot(sorted(u), p(u), label = 'p(u)', color = 'dodgerblue')
ax1.set_title('p(u)')

ax1.axis([0.98,1.02,p(u).min() - 8e-14,p(u).max() + 8e-14])

ax2.plot(sorted(u), q(u), label = 'q(u)', color = 'crimson')
ax2.set_title('q(u)')

ax2.axis([0.98,1.02,p(u).min() - 8e-14,p(u).max() + 8e-14])

ax1.set_xlabel('u')
ax2.set_xlabel('u')

plt.tight_layout()

plt.show()
