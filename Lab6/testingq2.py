#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from math import *

g = 9.81
l = 0.1

def f(r, t):
	theta = r[0]
	omega = r[1]
	ftheta = omega
	fomega = -(g/l)*sin(theta)

	return np.array([ftheta, fomega], float)

