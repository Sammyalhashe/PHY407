#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""
Using RK4, solves and plots the one dimensional 
system of the N identical masses with springs
"""

__author__ = "Eric Yeung"

m = 1.
k = 6.
omega = 2.
N = 5 # 5, 10, 20, 30
h = 0.001

# Define the driving forces of masses i

def F_i(i, t):

    if i == 0:
        drivingforce = np.cos(omega*t)

    else:
        drivingforce = 0

    return drivingforce

def f(r, t):
    
    # zeta1, zeta2, ..., zetaN, v1, v2, ..., vn = r
    
    zeta = r[:N] # Entrys before N 
    v = r[N:]    # Entrys after N

    # Let this be the first set of DEs
    dzeta = v

    # Create an empty list to store the second set of DEs with special cases
    dv = []
    
    for i in range(N):
        
        # Remember to shift indices, since python starts with 0
        # The starting point where the driving force F_i is non-zero
        if i == 0:
            dv.append(k/m*(zeta[1] - zeta[0]) + F_i(0, t)/m) 

        # The endpoint where the driving force F_N = 0
        elif i == N-1:
            dv.append(k/m*(zeta[N-2] - zeta[N-1]) + F_i(N, t)/m)

        # Every other value in between, with F_i = 0
        else:
            dv.append(k/m*(zeta[i+1] + zeta[i-1] - 2*zeta[i]) + F_i(i, t)/m)

    # Convert to arrays
    dzetaArray = np.array(dzeta, float)
    dvArray = np.array(dv, float)
    fArray = np.hstack((dzetaArray, dvArray)) # combine the two arrays horizontally
    
    return fArray 

times = np.arange(0, 20, h)

# Initial value for r

r = np.zeros(2*N)

zeta_sol = np.zeros((N, len(times)))

# Use runge-kutta method, and intialize t_index

def RK4(r, f):
    t_index = 0
    for t in times:
        zeta_sol[:, t_index] = r[:N] # Set each column to the r value  

        k1 = h*f(r, t)
        k2 = h*f(r + 0.5*k1, t + 0.5*h)
        k3 = h*f(r + 0.5*k2, t + 0.5*h)
        k4 = h*f(r + k3, t + h)

        r += (k1 + 2*k2 + 2*k3 + k4)/6
        t_index += 1
        
    return r

# Call the function

RK4(r, f) 

# Simple plot

for n in range(N):
    line = plt.plot(times, zeta_sol[n,:], label = 'M_' + str(n+1))    

plt.xlabel('time')
plt.ylabel('Amplitude')
plt.title(str(N) + 'Masses with Coupled Springs: $\\xi_i$(t)')

"""
# Animate

plt.ion()
for t in range(10):
    test = zeta_sol[n,:]
    line[0].set_ydata(test)
    plt.pause(0.01)

scat = plt.scatter(times, zeta_sol[:, 0] color = 'r')
line = plt.plot(times, zeta_sol[:, 0] )
"""

plt.legend().draggable()
plt.grid(True)
plt.show()
