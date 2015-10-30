#!/usr/bin/env python
from __future__ import division
from math import *

'''lab1_q1c.py: Determining and plotting the scattering angle distribution'''


__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

import matplotlib.pyplot as plt
import numpy as np

N = 2000 # number of particles
z = np.random.uniform(-1,1,N) # create a random array with N elements
t = map(lambda z:pi-2*asin(z), z) # apply our formula to the random values of height
tr = np.array(t) # make t into an array


# Count the values in our tr array for both ranges specified
FirstRange = ((170*pi/180 < tr) & (tr < 190*pi/180)).sum()
SecondRange = ((90*pi/180 < tr) & (tr < 110*pi/180)).sum()

print "%s scattering angles lie between 170 degrees and 190 degrees." % FirstRange
print "%s scattering angles lie between 90 degreesand 110 degrees." % SecondRange

relativeprob1 = FirstRange/N
relativeprob2 = SecondRange/N

print "The relative probability of finding the particle in the range(170,190) is %s for N = %s" % (relativeprob1, N)

print "The relative probability of finding the particle in the range(90,110) is %s for N = %s" % (relativeprob2, N)

import scipy.stats as stats

# Plot a gaussian to compare to uniform dist. 

gaussian = stats.norm.pdf(sorted(tr), np.mean(sorted(tr)), np.std(sorted(tr))) 
plt.plot(sorted(tr), gaussian, color = 'red', ls = '-', linewidth = '1.5')

# Plot the histogram with 100 bins

plt.hist(tr, bins = 100, color='b', normed = True)

plt.xlabel('Scattering Angle ($\Theta$)')
plt.title('Distribution of $\Theta$')

plt.show()