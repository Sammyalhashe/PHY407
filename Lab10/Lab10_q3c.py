#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import random 

"""
Compares 100 different iterations of importance sampling
and plots as a histogram
"""

__author__ = "Eric Yeung"

def f(x):
	return x**(-0.5)/(np.exp(x)+1) # Singularity at x = 0

def w(x):
	return x**(-0.5)

N = 10000
sum_fwx = 0.    
sum_fwxsq = 0.

a = 0.
b = 1.

# Store integration values in this array (100 values)
Iarray2 = np.zeros(100)

for j in range(100):
	for i in range(N):
		"""

		p(x) = 1/(2*sqrt(x)) -> int(p(x)) from 0 to x(z) = int(z) from 0 to z = z
		int(p(x)) = sqrt(x(z)) = z 
		---> x(z) = z**2
		Showed this in the writeup
		"""

		x = ((b-a)*random.random()+a)**2
		fwx = f(x)/w(x)

		sum_fwx += fwx
		sum_fwxsq += fwx**2

	# Divive sum_fx and sum_fxsq to for each (average) value in Iarray
	sum_fwx /= N
	sum_fwxsq /= N

	# From the lecture, I = <f(x)/w(x)> * integral of w(x) which is just 2 analytically
	Iarray2[j] = 2*sum_fwx

if __name__ == "__main__": 
	plt.hist(Iarray2,10,range=[0.8, 0.88])
	plt.xlabel('Importance Sampling Result')
	plt.title('Comparing different Iterations of the Integral')
	plt.show()
