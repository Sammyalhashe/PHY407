#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from Lab6_q2a import *

"""
Animates the vibrations of the masses
"""

__author__ = "Eric Yeung"

plt.ion()

for t in range(500):
    y = RK4(r,f)
    line[0].set_ydata(y)
    plt.title('Timestep = ' + str(t))
    plt.pause(0.001)

plt.show()