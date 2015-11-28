#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import random 

"""
Compares 100 different iterations of the mean value method
and plots as a histogram
"""

__author__ = "Eric Yeung"

def f(x):
	return x**(-0.5)/(np.exp(x)+1) # Singularity at x = 0

N = 10000
sum_fx = 0.    
sum_fxsq = 0.

a = 0.
b = 1.

# Store integration values in this array (100 values)
Iarray1 = np.zeros(100)

for j in range(100):
	for i in range(N):
		x = (b-a)*random.random()+a
		fx = f(x)

		sum_fx += fx
		sum_fxsq += fx**2

	# Divive sum_fx and sum_fxsq to for each (average) value in Iarray
	sum_fx /= N
	sum_fxsq /= N

	"""
	DO NOT USE THIS MEAN (Only first value is correct)
	mean = sum_fx/N  
	var = sum_fxsq/N - mean**2
	"""

	Iarray1[j] = (b-a)*sum_fx

if __name__ == "__main__": 
	plt.hist(Iarray1,10,range=[0.8, 0.88])
	plt.xlabel('Mean Value Method Result')
	plt.title('Comparing different Iterations of the Integral')
	plt.show()
