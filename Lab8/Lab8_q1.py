from __future__ import division, print_function
from pylab import *
# Constants
M = 100         # Grid squares on a side
V = 1.0         
target = 1e-6   # Target accuracy

# Create arrays to hold potential values
phi = zeros([M+1,M+1],float)

# Main loop
delta = 1.0

while delta>target:
	
	delta = 0
	for i in range(1,M):
		for j in range(1,M):
			if 20<i<80 and j==20:
				phi[i,j] = V
			elif 20<i<80 and j==80:
				phi[i,j] = -V
			else:
				difference =  (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1]) /4 - phi[i,j]
				phi[i,j] += difference
			if difference > delta:
				delta = difference
imshow(phi)
gray()
show()
