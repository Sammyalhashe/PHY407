#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""
Finding roots using various methods.
"""

__author__ = "Eric Yeung"

# Constants
Tmax = 2.0
points = 1000
accuracy = 1e-6

# Set up lists for plotting
y = []
temp = np.linspace(0.01,Tmax,points)

# Temperature loop
for T in temp:
    m1 = 1.0
    error = 1.0

    # Loop until error is small enough
    while error>accuracy:
        m1,m2 = np.tanh(m1/T),m1
        error = abs((m1-m2)/(1-T*np.cosh(m1/T)**2))
    y.append(m1)

# Make the graph
plt.plot(temp,y)
plt.ylim(-0.1,1.1)
plt.xlabel("Temperature")
plt.ylabel("Magnetization")

plt.show()

# Newtons method

yn = []

def f(m, T):
    return np.tanh(m/T) - m # rearrange

def deltam(m):
    return ((f(m + delta, T)) - f(m, T))/delta

def newton(f, m, T, delta = 1e-15):
    for T in temp:
        m1 = 1.
        error = 1.

        while error > accuracy:

            """
            if error < accuracy:
                break
            """

            m1,m2 = m - f(m1)/deltam(m1), m - f(m2)/deltam(m2)
            error = abs((m1-m2)/(1-T*np.cosh(m1/T)**2))
            print m1, m2
        yn.append(m1)

    return m1

print yn

"""
plt.plot(temp,yn)
plt.ylim(-0.1,1.1)
plt.xlabel("Temperature")
plt.ylabel("Magnetization")
plt.title('Newton\'s Method')

plt.show()

"""