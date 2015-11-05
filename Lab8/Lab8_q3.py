#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Calculates the temperature profile of the crust as a function 
of depth up to 20 metres and time up to 10 years.
Modified heat.py from textbook. 
"""

__author__ = "Eric Yeung"

# Constants
L = 20        # Depth in metres
D = 0.1       # Thermal diffusivity (m^2/day)
N = 100       # Number of divisions in grid
a = L/N       # Grid spacing
h = 1e-3     # Time-step
epsilon = h/100

Tmid = 10.0   # Temperature everywhere except surface and 20m
Thi = 11.0    # Temperature at 20m, deepest point

# Temperature at the surface
def T0(t):
	A = 10
	B = 12
	tau = 365
	return A + B*np.sin(2*np.pi*t/tau)

# Times to plot the temperature profile for the last year in three month intervals
# 10 years so multiply by 10, plot last year
"""
t1 = epsilon
t2 = t1 + 3376
t3 = t2 + (30 + 31 + 30)
t4 = t3 + (31 + 31 + 30)
t5 = t4 + (31 + 30 + 31)
tend = t5 + epsilon
"""
t1 = 365*9
t2 = t1 + (30 + 31 + 30)
t3 = t2 + (30 + 31 + 30)
t4 = t3 + (31 + 31 + 30)
t5 = t4 + (31 + 30 + 31)

# Create arrays
T = np.empty(N+1,float)
T[0] = Thi
T[N] = T0(epsilon)
T[1:N] = Tmid
Tp = np.empty(N+1,float)
Tp[0] = Thi
Tp[N] = T0(epsilon)

# Main loop
t = 0.0
c = h*D/(a*a)

while t < t1:

    Tp[1:N] = T[1:N] + c*(T[0:N-1] + T[2:N+1] - 2*T[1:N])
    T, Tp = Tp, T
    t += h

    # Make plots at the given times
    if abs(t - t1) < epsilon:
        plt.plot(T)
    if abs(t - t2) < epsilon:
        plt.plot(T)
    if abs(t - t3) < epsilon:
        plt.plot(T)
    if abs(t - t4) < epsilon:
        plt.plot(T)
    if abs(t - t5) < epsilon:
        plt.plot(T)

plt.xlabel("Depth")
plt.ylabel("Temperature")
plt.title("Thermal Diffusion in Earth's Crust")

plt.show()
