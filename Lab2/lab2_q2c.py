#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from scipy import special as sp

"""
Computes the integral for bessel function of first kind, and compares
to the bessel function in scipy. The two are plotted against each other.
Then makes a density plot. 
"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

# m is nonnegative integer, x >= 0
N = 1000

def J(m,x):
    
    def f(theta):
        return cos(m*theta - x*sin(theta))

    a = 0
    b = pi
    h = (b-a)/N

    s = f(a) + f(b) # from 5.9
    for i in range(1, N, 2):  # odd terms
        s += 4*f(a+i*h)
    for i in range(2, N-1, 2): # even terms
        s += 2*f(a+i*h)
    
    return (1/pi)*(h/3)*s

print J(0,1), sp.jv(0,1)  # Intial comparison test


x = np.linspace(0,20)

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True)

ax1.plot(x, J(0,x), color = 'r', label = 'J(0,x)')
ax1.plot(x, J(1,x), color = 'b', label = 'J(1,x)')
ax1.plot(x, J(2,x), color = 'g', label = 'J(2,x)')

ax1.set_xlabel('x')
ax1.set_ylabel('J(m,x)')
ax1.set_title('$J_0$, $J_1$, and $J_2$: Homemade')
ax1.legend(loc = 'upper right')


ax2.plot(x, sp.jv(0,x), color = 'r', label = 'J(0,x)')
ax2.plot(x, sp.jv(1,x), color = 'b', label = 'J(1,x)')
ax2.plot(x, sp.jv(2,x), color = 'g', label = 'J(2,x)')

ax2.set_xlabel('x')
ax2.set_ylabel('J(m,x)')
ax2.set_title('$J_0$, $J_1$, and $J_2$: Scipy')
ax2.legend(loc = 'upper right')

plt.tight_layout()
plt.show()


###########################
######  Density Plot ######
###########################

# Define intensity as a function of x and y
# Referred to page 109 in the textbook for this part. 

def I(x, y):
    r = sqrt(x**2+y**2)
    if (r == 0):
        return 1/2   # Using the limit found in the textbook. 
    else:
        return (sp.jv(0,k*r)/(k*r))**2

wavelength = 500   # in nano meters
k = 2*pi/wavelength # our wavenumber 

# Value in which we increase our x, and y
step = 0.95 

# Create an empty matrix
Imatrix = np.zeros(shape=(N,N))

# Calculate the values in the matrix above
for i in range(N):
    y = step*i - 500 # 500nm shifted from top, should center
    for j in range(N):
        x = step*j - 500 # 500 nm shifted from right side, should center
        Imatrix[i,j] = I(x,y)
            
plt.imshow(Imatrix, vmax= 0.04, cmap='hot')
plt.xlabel('x (nm)')
plt.ylabel('y (nm)')
plt.title('Diffraction Pattern in Telescopes')

plt.colorbar() 
plt.show()