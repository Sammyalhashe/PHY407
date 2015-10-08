#!/usr/bin/env python
from __future__ import division
from math import *
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

"""
Displays the roundoff errors in a product f(u) instead of series
"""

__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

N = 100
u = np.random.uniform(0.98,1.02, N)

f = u**8/((u**4)*(u**4))
newf = f-1


plt.scatter(u, newf, color = 'r')

#print newf.min(), newf.max()
plt.ylim([newf.min()+1e-16, newf.max()+1e-17]) # Change axis to view more dots
plt.xlim([u.min(), u.max()]) # Same as above

plt.xlabel('u')
plt.ylabel('f(u)-1')

plt.title('The roundoff errors in f-1')
plt.show()

print np.mean(f), np.std(f)

C = 1e-16

estimateError = sqrt(2)*C*1**8/((1**4)*(1**4)); print estimateError