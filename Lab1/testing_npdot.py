#!/usr/bin/env python
from __future__ import division
from math import *
import numpy as np

import time

for Z in np.arange(2, 301): # 2-300 calculations, slow(nested loops, vectorize?)

	A = np.ones([Z,Z],float)*3.
	B = np.ones([Z,Z],float)*3.
	C = np.zeros([Z,Z],float)

	# Now let's time the np.dot function

	startTime2 = time.time()

	D = np.dot(A, B)

	time2 = time.time() - startTime2
	print "Time for np.dot for %sx%s matrix is %s" % (Z, Z, time2)
