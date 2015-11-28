#!/usr/bin/env python
from __future__ import division
import numpy as np
import random 

"""
Computes integral using mean value monte carlo method
with 1e4 points
"""

__author__ = "Eric Yeung"

def f(x):
	return x**(-0.5)/(np.exp(x)+1) # Singularity at x = 0

N = 10000
sum_fx = 0.    
sum_fxsq = 0.

a = 0.
b = 1.

for i in range(N):
	x = (b-a)*random.random()+a
	fx = f(x)

	sum_fx += fx
	sum_fxsq += fx**2

mean = sum_fx/N
var = sum_fxsq/N - mean**2

I = (b-a)*mean
sigma = (b-a)*np.sqrt(var/N)

print "Mean value method gives I = %s" % I
