#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

def p(u):
    return (1-u)**8

def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

N = 1000 # number of points
u = np.random.uniform(0.98,1.02,N)

import scipy.stats as stats

# Plot a gaussian to compare to uniform dist. 

gaussian = stats.norm.pdf(sorted(p(u)-q(u)), np.mean(sorted(p(u)-q(u))), np.std(sorted(p(u)-q(u)))) 
plt.plot(sorted(p(u)-q(u)), gaussian, color = 'red', ls = '-', linewidth = '1.5', label ='Gaussian')

plt.hist(p(u)-q(u), bins = 100, color='g', normed = True, label ='p(u)-q(u)')
plt.xlabel('p(u)-q(u)')
plt.title('p(u)-q(u) distribution: u near 1')

plt.legend(loc='upper left')

plt.show()

###############################################################
# Calculating the standard deviation manually, and with numpy # 
###############################################################

# Create an array of terms

qArray = [1, -8*u, 28*u**2,- 56*u**3, 70*u**4, -56*u**5, 28*u**6, -8*u**7, u**8]

C = 1e-16

# For loop that squares each term of qArray

for k in range(9):
    qsquared = qArray[k]**2
    k += 1
    
# Calculate the mean using generic formula

calcmean = 1/(len(qsquared)) * np.sum(qsquared)

# Now calculate sigma

calcsigma = C * sqrt(len(qsquared))*sqrt(calcmean) 

# Compare

print calcsigma, np.std(p(u)-q(u))