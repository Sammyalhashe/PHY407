#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

"""
Density plot of values of the photo
"""

__author__ = "Eric Yeung"

bArray = np.loadtxt("blur.txt"); # print bArray.shape
plt.imshow(bArray, cmap = 'gray')
plt.show()