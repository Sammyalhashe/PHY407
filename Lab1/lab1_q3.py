#!/usr/bin/env python
from __future__ import division
from math import *

'''lab1_q3.py: Demonstrating the relationship between matrix size and runtime'''


__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

import matplotlib.pyplot as plt
import numpy as np
import time


for Z in np.arange(2, 301): # 2-300 calculations, slow(nested loops, vectorize?)

	# Create the matrices A, B full of 3s
	# Create an empty matrix C

	A = np.ones([Z,Z],float)*3.
	B = np.ones([Z,Z],float)*3.
	C = np.zeros([Z,Z],float)

	# Start timing first operation

	startTime1 = time.time()

	for i in range(Z):
	    for j in range(Z):
	        for k in range(Z):
	            C[i,j] += A[i,k]*B[k,j]


	time1 = time.time() - startTime1 # Stop timing first operation 

	print "Time for normal matrix multiplication of %sx%s matrix is %s" % (Z, Z, time1)

	# Now let's time the np.dot function

	startTime2 = time.time()

	D = np.dot(A, B)

	time2 = time.time() - startTime2 # Stop timing second operation
	
	print "Time for np.dot for %sx%s matrix is %s" % (Z, Z, time2)

	plt.scatter(Z, time1, color = 'b', label = "norm")
	plt.scatter(Z, time2, color = 'r', label = "npdot")

	plt.xlabel('Z')
	plt.ylabel('Runtime')
	plt.title('Relationship between Z[2,300] and runtime')

#plt.legend(loc='upper right') don't need a legend