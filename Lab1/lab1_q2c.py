#!/usr/bin/env python
from __future__ import division
from math import *

'''lab1_q2c.py: Change of intial conditions: x_c'''


__author__  = "Eric Yeung"
__email__   = "eric.yeung@mail.utoronto.ca"


import matplotlib.pyplot as plt
import numpy as np
from pylab import *

# Now we repeat the previous question! 

# New parameters
xc = 8.654e7 # as shown in the pdf explanation  (in meters)
dt = 1.0e-5
# Old parameters
k = 12.
m = 1. # in metres 
c = 2.998e8
ti = 0.
tf = 10.

time=arange(ti,tf,dt)
tlen=len(time)

# First we do the example using x_c/10 as the initial displacement

x = np.zeros(tlen)
v = np.zeros(tlen)

# initial conditions

x[0] = xc
v[0] = 0.

i = 0

# Classical particle, so no relativistic effects, use that system of ODEs
while i < tlen-1: 
    x[i+1] = x[i] + v[i]*dt
    v[i+1] = v[i] - k/m*x[i]*dt
    i += 1


# Now onto the example using x_c as the initial displacement 

xnew = np.zeros(tlen)
vnew = np.zeros(tlen)    
    
# Different initial conditions 

xnew[0] = 10*xc
vnew[0] = 0


# Repeat but instead of i, use different index j 
j = 0

while j < tlen-1: 
    xnew[j+1] = xnew[j] + vnew[j]*dt
    vnew[j+1] = vnew[j] - k/m*xnew[j]*dt
    j += 1
    
f, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, sharex='row', sharey='col')
ax1.plot(time, x, color='springgreen')
ax1.set_title('Position: xc')

ax2.plot(time, xnew, color='darkgreen')
ax2.set_title('Position: 10*xc')

ax3.plot(time, v, color='gold')
ax3.set_title('Velocity: xc')

ax4.plot(time, vnew, color='darkgoldenrod')
ax4.set_title('Velocity: 10*xc')

ax1.set_xlabel('Time')
ax2.set_xlabel('Time')
ax3.set_xlabel('Time')
ax4.set_xlabel('Time')

ax1.set_ylabel('Position')
ax2.set_ylabel('Position')
ax3.set_ylabel('Velocity')
ax4.set_ylabel('Velocity')

plt.tight_layout()
#plt.savefig('new_graphs.png', dpi=4000, bbox_inches="tight")
plt.show()

print x[999999] - xnew[999999] # Let the system evolve until the limit (1e6)
print v[999999] - vnew[999999] 
# Pretty big differences here. The only changes I see are the amplitudes. 