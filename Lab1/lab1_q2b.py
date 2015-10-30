#!/usr/bin/env python
from __future__ import division
from math import *

'''lab1_q2b.py: Comparing the relativistic and non-relativistic mass on a spring'''


__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"

import matplotlib.pyplot as plt
import numpy as np
from pylab import *

# Define some parameters

k = 12.
m = 1. 
c = 2.998e8
dt = 0.001  # timestep
ti = 0.
tf = 10.

time=arange(ti,tf,dt) # Define time array
tlen=len(time) # Gives us the dimension/length of the array

# Create empty arrays of dimension tlen 

x = zeros(tlen)
v = zeros(tlen)

# initial conditions

x[0] = 1
v[0] = 0

# Initial index value

i = 0

# Euler Method

while i < tlen-1:
    x[i+1] = x[i] + v[i]*dt
    v[i+1] = v[i] - k/m*x[i]*(1- (v[i])**2/c**2)**(3/2)*dt
    i += 1


# With no relativistic effects, the v^2/c^2 -> 0 in our ODE. This makes our system a lot cleaner! 
# Define new variables xn, vn for non-relativistic position and velocities respectively. 

xn = zeros(tlen)
vn = zeros(tlen)

xn[0] = 1
vn[0] = 0

# Different indices, so no confusion

j = 0 

# Euler Method again, straightforward

while j < tlen-1: 
    xn[j+1] = xn[j] + vn[j]*dt
    vn[j+1] = vn[j] - k/m*xn[j]*dt
    j += 1


# Show explicitly the difference after time evolution by comparing last 
# entrys in the arrays x, xn, v, and vn 

print x[9999] - xn[9999] # Let the system evolve until the limit (1e4)
print v[9999] - vn[9999]

# Plot everything

f, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, sharex='row', sharey='col')
ax1.plot(time, x, color='deepskyblue')
ax1.set_title('Position Graph: Relativistic')

ax2.plot(time, xn, color='navy')
ax2.set_title('Position Graph: Non-Relativistic')

ax3.plot(time, v, color='crimson')
ax3.set_title('Velocity Graph: Relativistic')

ax4.plot(time, vn, color='darkred')
ax4.set_title('Velocity Graph: Non-Relativistic')

ax1.set_xlabel('Time')
ax2.set_xlabel('Time')
ax3.set_xlabel('Time')
ax4.set_xlabel('Time')

ax1.set_ylabel('Position')
ax2.set_ylabel('Position')
ax3.set_ylabel('Velocity')
ax4.set_ylabel('Velocity')

plt.tight_layout()
#plt.savefig('NR_R_graphs.png', dpi=4000, bbox_inches="tight")
plt.show()
