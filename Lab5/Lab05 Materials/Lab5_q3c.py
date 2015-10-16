#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

"""
Unblurs the photo generated from blur.txt
"""

__author__ = "Eric Yeung"

sigma = 25
epsilon = 1e-3

# Read data from blur.txt

bArray = np.loadtxt("blur.txt")
rows = bArray.shape[0]
cols = bArray.shape[1]

# Calculates the point spread function

def Gauss(sigma):
	
	GaussArray = np.empty([rows,cols], float)

	for i in range(rows):	
		ip = i
		if ip > rows/2:
			ip -= rows #bottom half of rows moved to negative values		
		for j in range(cols):
			jp = j
			if jp > cols/2:
				jp -= cols #right half of columns moved to negative values	
			GaussArray[i,j] = np.exp(-(ip**2+jp**2)/(2.0*sigma**2)) #compute gaussian

	return GaussArray

# Fourier Transform both distributions

Blurry = np.fft.rfft2(bArray)
fGauss = np.fft.rfft2(Gauss(sigma))

def division(epsilon):

	# Pick out the entrys that are too small for division

	smallentrys = np.where(abs(fGauss) < epsilon)
	bigentrys = np.where(abs(fGauss) >= epsilon)

	# Create an empty matrix to store divisions

	div = np.empty([rows, len(Blurry[0])], complex)

	div[bigentrys] = Blurry[bigentrys]/(fGauss[bigentrys]) # divide by all the big entrys 			
	div[smallentrys] = Blurry[smallentrys] # don't divide by small entrys or divide by "1"

	return div

# Take the inverse fourier transform and divide by M and N

NotBlurry = np.fft.irfft2(division(epsilon)/(rows*len(Blurry[0])))

# Plot and compare

if __name__ == '__main__':
	
	plt.imshow(Gauss(sigma), vmax = 0.005, cmap = 'gray')
	
	f, (ax1, ax2) = plt.subplots(1, 2)

	ax1.imshow(bArray, cmap = 'gray')
	ax1.set_title('Blurry image')

	ax2.imshow(NotBlurry, cmap ='gray')
	ax2.set_title('Unblurry image')

	plt.show()
