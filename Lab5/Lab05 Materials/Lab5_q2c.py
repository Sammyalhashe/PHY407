#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

"""
Calculates the velocity
"""

__author__ = "Eric Yeung"

SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

# From the formula in the lab
circumference = 2*np.pi*6.371e6*np.cos(50*np.pi/180); #print circumference 

dlatidude = Longitude*np.cos(np.pi*50/180); #print dlatidude

seconds = Times*86400; #print seconds

propagationspeed = (dlatidude[-1] - dlatidude[0])/(seconds[-1] - seconds[0]); print propagationspeed