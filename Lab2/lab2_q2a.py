#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy

"""
Comparing simpson's and trapezoidal rule for a generic integral
"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"


def f(x):
    return x**4 - 2*x + 1

# TRAPEZOIDAL 

N = 10
a = 0.0
b = 2.0
h = (b-a)/N

s1 = 0.5*f(a) + 0.5*f(b)

for k in range(1,N):
    s1 += f(a+k*h)

I1 = s1*h

print I1

# SIMPSONS

N = 10
a = 0.0
b = 2.0
h = (b-a)/N

s2 = f(a) + f(b)
    
for k in range(1,N,2):
    s2 += 4*f(a+k*h)

for k in range(2,N-1,2):
    s2 += 2*f(a+k*h)
 
I2 = s2*h/3

print I2


# Fractional Error
expectedvalue = 4.4
C = 1e-16

fError = abs(I2 - expectedvalue)/abs(expectedvalue)

print fError

