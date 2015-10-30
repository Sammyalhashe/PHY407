#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

"""
Solving comit orbits using RK4 with fixed step size

Our two equations become 4 equations:
	dx/dt = v_x
	dy/dt = v_y
	dv_x/dt = -G*M*x/radius**3
	dv_y/dt = -G*M*y/radius**3
"""

__author__ = "Eric Yeung"

M = 1 # Let this be in solar masses for convenience
G = 39.5 # In Au^3 * solarmass^-1 * yr^-2 (6.67408e-11 in m^3 * kg^-1 * s^-2)

h = 5e-4 # Choose our time step so the orbits lie on top of each other
int_time = 200 # 200 years

startTime1 = time.time()

# There is no dependence on t anywhere here, remove it as an argument
def f(r):
	x = r[0]
	y = r[1]
	v_x = r[2]
	v_y = r[3]
	radius = np.sqrt(x**2 + y**2)
	
	# First set of ODEs
	dx = v_x
	dy = v_y

	# Second set of ODEs
	dv_x = -G*M*x/radius**3
	dv_y = -G*M*y/radius**3
	
	fArray = np.array([dx, dy, dv_x, dv_y], float)
	return fArray

def RK4(r, f):
	k1 = h*f(r)
	k2 = h*f(r + 0.5*k1)
	k3 = h*f(r + 0.5*k2)
	k4 = h*f(r + k3)

	r += (k1 + 2*k2 + 2*k3 + k4)/6
	return r

# Define the intial R vector with our initial conditions

x0 = 26.74 # In AU (4e12 metres)
y0 = 0. # Somewhere near Neptune's orbit
vx0 = 0. 
vy0 = 0.1054 # In AU/yr (500 m/s)

r0 = np.array([x0, y0, vx0, vy0],float)

# Initialize an array to store our solutions

r_sols = np.empty([4, int_time/h])

# Initialize loop index
t = 0

while t*h < int_time: 
	rtemp = RK4(r0, f)
	r_sols[:, t] = rtemp
	# Set r0 as the new value at the new step
	r0 = rtemp
	t += 1

print r_sols[0,:], r_sols[1,:] # Debugging

# Plot x vs y
plt.plot(r_sols[0,:], r_sols[1,:])
plt.annotate('SUN', xy=(-1,0), bbox=dict(boxstyle="circle", fc="r"))
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comet Orbits')
plt.show()

time1 = time.time() - startTime1 # Stop timing operation 

if __name__ == "__main__": 
	print "Comet orbit solution took %s seconds" % time1
	print "The stepsize used was %s" % h

else:
	print "Importing from the fixed RK4..."