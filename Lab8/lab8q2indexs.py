#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Simulates linear water waves
"""

__author__ = "Eric Yeung"

L      = 400. # Metres
D      = 50. # Metres 
bdelta = 2. # Metres
A      = 1. # Metre
M      = 100
g 	   = 9.8 # m^2/s
target = 1e-4

xpoints = np.arange(0, L, L/M)
#zpoints = np.arange(-D, 0, D/M)
zpoints = np.arange(0, D, D/M)

eta = np.zeros(M, float)

# Gaussian so it decays exponentially away from the peak
for i in range(len(xpoints)):
    eta[i] = -A*np.exp(-(xpoints[i]-L/2)**2/bdelta**2)

# Initialize phi
phi = np.zeros([M, M], float)

"""
Euler Forward
"""

dx = 1.0
dz = 1.0
h  = 0.1

etafull = np.copy(eta)
etahalf = np.zeros(M, float)

phifull = np.copy(phi)
phihalf = np.zeros([M, M], float)

Initialdphi = (phifull[0,:] - phifull[1,:])/dz # dphi/dz at t = 0
phihalf[0,:] = phifull[0,:] - h/2*g*etafull

etahalf = etafull + h/2*Initialdphi # Second term should be zero

"""
Jacobi Method to solve Laplace's equation
"""

phihalfprime = np.empty([M, M], float)
delta = 1.

while delta > target:
	phihalfprime = np.copy(phihalf)
	phihalfprime[1:-1,1:-1]=(phihalfprime[2:,1:-1]+phihalfprime[:-2,1:-1]
	                       +phihalfprime[1:-1,2:]+phihalfprime[1:-1,:-2])/4

	delta = np.amax(abs(phihalf-phihalfprime))
	phihalf,phihalfprime = phihalfprime,phihalf

print phihalf

halfdphi = (phihalf[0,:] - phihalf[1,:])/dz # dphi/dz at t = h/2
nsteps = 200
i = 1

"""
while i < nsteps:
	# From h/2 to the next full step h
	etahalf += h*dx
	phihalf[0,:] += h*dz
	t +=  0.5*h
"""


"""
zpoints = np.arange(-D, 0, D/M)
Initialdphi = (phifull[-1,:] - phifull[-2,:])/dz # dphi/dz at t = 0

phihalf[-1,:] = phifull[-1,:] - h/2*g*etafull
"""
