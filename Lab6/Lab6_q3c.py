#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

"""
Solving comit orbits using RK4 with adaptive step size

	From page 359 of the textbook, the errors are:
	xerror = (x1 - x2)/30.
	yerror = (y1 - y2)/30.

and rho, from page 358, is:
	rho = 30*h*delta/abs(xerror + yerror)
"""

__author__ = "Eric Yeung"

M = 1 # Let this be in solar masses for convenience
G = 39.5 # In Au^3 * solarmass^-1 * yr^-2 (6.67408e-11 in m^3 * kg^-1 * s^-2)

h0 = 5e-2 # This is now our intial stepsize (5e-4)

int_time = 50 # 200 years

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

# Change the runge kutta function to a function of h.

def RK4(h, r, f):
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

r_sols = np.empty([4, int_time/h0])

# Set a target accyracy of delta = 1 km/yr

delta = 6.685e-9 # In Au/yr



"""
______________DEBUGING__________________________________

rtemp = RK4(h, r0, f)

# Error at h and 2*h
rError1 = RK4(h, rtemp, f)
rError2 = RK4(2*h, r0, f)

xerror = (rError1[0] - rError2[0])/30. # This is 0
yerror = (rError1[1] - rError2[1])/30. # This is 0 too

# This means rho is NAN always

print rError1; print rError2; print xerror; print yerror
________________________________________________________
"""

# Initialize loop index
t = 0
h = h0

while t < int_time: 
	rtemp = RK4(h, r0, f)
	
	# Error at h and 2*h
	rError1 = RK4(h, rtemp, f)
	rError2 = RK4(2*h, r0, f)
	xerror = (rError1[0] - rError2[0])/30.  # 0th element -> x component
	yerror = (rError1[1] - rError2[1])/30.
	rho = h*delta/np.sqrt(xerror**2 + yerror**2)

	# If rho > 1, actual accuracy is better than the target accuracy. Keep it.
	if rho > 1:
		t += h
		r0 = rtemp
		#h = h*rho**(1/4) # Make it bigger since rho^1/4 > 1
		r_sols[:, t] = rtemp	

	# If rho < 1, not small enough, h' = h*rho^(1/4) from newman, discard it.
	elif rho < 1:
		h = h * rho**(1/4)

	# If denominator is very small, let us just default back to our initial stepsize 
	#elif abs(rError1[0] - rError2[0]) < 1e-14:
	#	h = h0 

	# Any other case not listed
	else:
		h = 10*h0 # Need an upper bound

time1 = time.time() - startTime1 # Stop timing operation 

if __name__ == "__main__": 
	# Plotting x vs y
	plt.plot(r_sols[0,:], r_sols[1,:], '-.')
	plt.annotate('SUN', xy=(-1,0), bbox=dict(boxstyle="circle", fc="r"))
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Comet Orbits: Adaptive Step Size')
	plt.show()

	print "Adaptive Comet orbit solution took %s seconds" % time1
	
else:
	print "Importing from the adaptive step size RK4 program..."
