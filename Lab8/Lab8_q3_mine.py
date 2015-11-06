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
h = 1e-2      # Time-step
epsilon = h/100

# Temperature at the surface
def T0(t):
    A = 10
    B = 12
    tau = 365  # 365 days, or 12 months
    return A + B*np.sin(2*np.pi*t/tau)

# Times to plot the temperature profile for the last year in three month intervals
# 10 years so multiply by 10, plot last year


"""

t1 = epsilon
t2 = t1 + 3376            # Start of 10th year
t3 = t2 + (30 + 31 + 30)
t4 = t3 + (31 + 31 + 30)
t5 = t4 + (31 + 30 + 31)
tend = t5 + epsilon       # End of 10th year

"""


t1 = epsilon
t2 = t1 + (30 + 31 + 28)   # Start of 10th year
t3 = t2 + (30 + 31 + 30)
t4 = t3 + (31 + 31 + 30)
t5 = t4 + (31 + 30 + 31)
tend = t5 + epsilon       # End of 10th year

timerange = np.arange(0, tend, tend/(N+1))
depth = np.arange(0, 20, 20/(N+1))

# Our intial conditions are now 2D, tmid and tdeep are BC same in entire row
Tsurface = np.zeros(len(timerange), float) 
for t in range(len(timerange)):
    Tsurface[t] = T0(t)

Tmid = np.zeros(len(timerange), float)
for t in range(len(timerange)):
    Tmid[t] = 10.0 # Temperature everywhere except surface and 20m INTIAL CONDITION ONLY

Tdeep = np.zeros(len(timerange), float)
for t in range(len(timerange)):
    Tdeep[t] = 11.0 # Highest temperature on the surface

# Create arrays
T = np.empty([N+1, N+1], float)
Tp = np.empty([N+1, N+1], float)

T[0,:] = Tsurface # Start with the temp at the surface, -2 lowest
T[N,:] = Tdeep  # Start with the temp at the deepest point, 11
T[1:N,:] = Tmid # Start at 10 degrees everywhere except endpoints
Tp[0,:] = Tsurface 
Tp[N,:] = Tdeep 

# Main loop
t = 0.0
c = h*D/(a*a)

while t < tend:

    #Tp[1:N,:] = T[1:N,:] + c*(T[0:N-1,:] + T[2:N+1,:] - 2*T[1:N,:])
    Tp[:,1:N] = T[:,1:N] + c*(T[:,0:N-1] + T[:,2:N+1] - 2*T[:,1:N])

    T, Tp = Tp, T
    t += h

    # Make plots at the given times
    if abs(t - t1) < epsilon:
        plt.plot(T[1])
    if abs(t - t2) < epsilon:
        plt.plot(T[1])
    if abs(t - t3) < epsilon:
        plt.plot(T[1])
    if abs(t - t4) < epsilon:
        plt.plot(T[1])
    if abs(t - t5) < epsilon:
        plt.plot(T[1])

plt.xlabel("Depth")
plt.ylabel("Temperature")
plt.title("Thermal Diffusion in Earth's Crust")

plt.show()

