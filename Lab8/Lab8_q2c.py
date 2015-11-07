#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Change L = 1500, nsteps = 300, and h = 0.5
"""

__author__ = "Eric Yeung"

L      = 1500.# Metres
D      = 50. # Metres 
bdelta = 2.  # Metres
A      = 1.  # Metre
g      = 9.8 # m^2/s
target = 1e-4 # Initial thershold

# This is our spacing/resolution
dx = 1.0
dz = 1.0

xpoints = np.arange(0, L, dx)
zpoints = np.arange(-D, 0, dz)

# Initialize eta (1D Array)
eta = np.zeros(len(xpoints), float)

# Initialize phi (2D Array)
phi = np.zeros([len(zpoints), len(xpoints)], float)

# Start timing 
import time
startTime1 = time.time()

# Gaussian so it decays exponentially away from the peak
for i in range(len(xpoints)):
    eta[i] = -A*np.exp(-(xpoints[i]-L/2)**2/bdelta**2)

"""
Euler Forward method to get our half step
"""

h = 0.5

etafull = np.copy(eta)
etahalf = np.zeros(len(xpoints), float)

phifull = np.copy(phi)
phihalf = np.zeros([len(zpoints), len(xpoints)], float)

# Finite differences dphi/dz at t = 0
fulldphi = (phifull[0,:] - phifull[1,:])/dz 

etahalf = etafull + h/2*fulldphi # Second term should be zero
phihalf[0,:] = phifull[0,:] - h/2*g*etafull

"""
Jacobi Method to solve Laplace's equation
"""

delta = 1.

def Laplace(delta, arg):
	while delta > target:
		phiprime = np.copy(arg)
		phiprime[1:-1,1:-1] = (phiprime[2:,1:-1] + phiprime[:-2,1:-1]
		                       +phiprime[1:-1,2:] + phiprime[1:-1,:-2])/4

		delta = np.amax(abs(arg - phiprime))
		arg, phiprime = phiprime, arg

	return arg

# Call the function to get phi(x, z=0, h/2)
phihalf = Laplace(delta, phihalf)

# Finite Differences for phi at t = h/2
halfdphi = (phihalf[0,:] - phihalf[1,:])/dz # dphi/dz at t = h/2

"""
Main Leapfrogging Loop
"""

nsteps = 300
i = 1
t = 0

# Saving the values for plotting
etavalues = []
etavalues.append(eta) # append initial conditions
phivalues = []

while i < nsteps:

	fulldphi = (phifull[0,:] - phifull[1,:])/dz
	halfdphi = (phihalf[0,:] - phihalf[1,:])/dz

	# From h/2 to the next full step h using Derivative at the half step
	etafull = etafull + h*halfdphi 
	phifull[0,:] = phifull[0,:] -h*g*etahalf

	# Solve Laplace's equation
	phifull = Laplace(delta, phifull)	

	t += h 

	fulldphi = (phifull[0,:] - phifull[1,:])/dz
	halfdphi = (phihalf[0,:] - phihalf[1,:])/dz

	# From h to 3h/2
	etahalf = etahalf + h*fulldphi
	phihalf[0,:] = phihalf[0,:] -h*g*etafull

	# Solve Laplace's equation again
	phihalf = Laplace(delta, phihalf)

	# Save eta and phi from each full step
	etavalues.append(etafull)
	phivalues.append(phifull)

	# Increment
	t += 0.5*h
	i += 1

	# Find out which index gives us the times we want (0, 2, 5)
	#print t, i

"""
Pick out the rows that correspond to the different columns 
and start to plot x vs eta 
"""

etaplot0 = np.array(etavalues)[0,:] # at t = 0   (i = 0)
etaplot2 = np.array(etavalues)[15,:] # at t = 2  (i = 41)
etaplot5 = np.array(etavalues)[35,:] # at t = 5 (i = 101)
etaplot10 = np.array(etavalues)[70,:] # at t = 5 (i = 101)
etaplot20 = np.array(etavalues)[135,:] # at t = 5 (i = 101)
etaplot30 = np.array(etavalues)[199,:] # at t = 5 (i = 101)

phiplot = np.array(phivalues)[0][0,:] # Don't need to plot phi anyways

figure, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex = True)
ax1.plot(xpoints, etaplot0, color = 'maroon', label = 't = 0s')
ax2.plot(xpoints, etaplot2, color = 'deepskyblue', label = 't = 2s')
ax3.plot(xpoints, etaplot5, color = 'forestgreen',label = 't = 5s')
ax4.plot(xpoints, etaplot10, color = 'darkorange',label = 't = 10s')
ax5.plot(xpoints, etaplot20, color = 'navy',label = 't = 20s')
ax6.plot(xpoints, etaplot30, color = 'fuchsia',label = 't = 30s')

time1 = time.time() - startTime1 # Stop timing second operation
print "This simulation took %s seconds!" % time1

# Hide y-tick labels because too messy
for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
	ax.set_yticklabels([])

figure.text(0.0, 0.5, '$\eta$(x,t)', va='center', rotation='vertical')
ax6.set_xlabel('x')

ax1.set_title("Linear Water Waves at t = (0, 2, 5, 10, 20, 30) s")
plt.tight_layout()

if __name__ == "__main__": 
	plt.show()