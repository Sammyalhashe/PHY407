#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""
Calculates the phonon dispersion for Pyrite, 
a lattice with a 3 atom basis
"""

__author__ = "Eric Yeung"

k = 5   # spring constant in N/m
m = 1   # mass in kg
a = 1.5   # interplanar spacing

kappa = np.linspace(-np.pi/a, np.pi/a)

omega = np.sqrt(4*k/m)*abs(np.sin(kappa*a/2))

plt.plot(kappa, omega)
plt.show()