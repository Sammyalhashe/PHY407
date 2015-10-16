#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

"""
Density plot of point spread function
"""

__author__ = "Eric Yeung"

bArray = np.loadtxt("blur.txt")

rows = bArray.shape[0]
cols = bArray.shape[1]

def Gauss(rows, cols, sigma):
	
	GaussArray = np.empty([rows,cols], float)

	for i in range(rows):	
		ip=i
		if ip>rows/2:
			ip -=rows #bottom half of rows moved to negative values		
		for j in range(cols):
			jp=j
			if jp>cols/2:
				jp -= cols #right half of columns moved to negative values	
			GaussArray[i,j]=np.exp(-(ip**2+jp**2)/(2.0*sigma**2)) #compute gaussian

	return GaussArray

plt.imshow(Gauss(rows, cols, sigma), vmax = 0.005, cmap = 'gray')
plt.show()