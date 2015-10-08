#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy

"""Calculates a generic integral using trapezoidal rule with N=10 and N=20"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

def f(x):
    return x**4 - 2*x + 1

N1 = 10
a = 0.0
b = 2.0
h1 = (b-a)/N1

s1 = 0.5*f(a) + 0.5*f(b)
for k in range(1,N1):
    s1 += f(a+k*h1)

print h1*s1

#############################

N2 = 20

h2 = (b-a)/N2

s2 = 0.5*f(a) + 0.5*f(b)
for j in range(1,N2):
    s2 += f(a+j*h2)

print h2*s2
   
ErrorEst = abs(1/3*(h2*s2 - h1*s1))

print ErrorEst

# Compare this error to the error between the last calculation and the true value of 4.4

Error1 = abs(h1*s1 - 4.4)/abs(4.4); print Error1

Error2 = abs(h2*s2 - 4.4)/abs(4.4); print Error2
