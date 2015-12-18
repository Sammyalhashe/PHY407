#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Variational Monte Carlo Method for the Helium Atom. 
We have 2 electrons, so 2 independent positions.
"""

__author__ = "Eric Yeung"

np.random.seed(1337) # Seed the random generator 

N = 1000 # Number of Monte Carlo steps
h = 1.00 # Step-size

# Trial wave function
def psi(alpha, r1, r2):
	r12 = abs(r1-r2) # Relative vector 
	return np.exp(-2*r1)*np.exp(-2*r2)*np.exp(r12/(2*(1+alpha*r12)))

# Local energy 
def localenergy(alpha, r1, r2):
	r12 = abs(r1-r2)
	return -4 + alpha/(1+alpha*r12) + alpha/(1+alpha*r12)**2 + alpha/(1+alpha*r12)**3 \
	- 1/(1+alpha*r12)**4 + r12*(r1-r2)/(1+alpha*r12)**2  

# Metropolis Algorithm to keep position or accept a trial position 
def metropolis(alpha, r1, r2, r1trial, r2trial):
	# Ratio between the probability functions |psi|^2
	omega = (psi(alpha, r1trial, r2trial)/psi(alpha, r1, r2))**2 

	if omega >= 1:
		return True # New position is accepted unconditionally
	else:
		if omega > np.random.rand(): # Acceptance probability
			return True # New position is accepted 
		else:
			return False # New position is not accepted, keep old 

# The main program
def main(alpha):
	# Radii, random numbers from 0 to 1
	r1 = np.random.rand()
	r2 = np.random.rand()

	# Initialize energy and energy squared
	E1 = 0.
	E2 = 0.

	# Keep track of steps
	steps = 0.

	# Begin the the loop for the Markov Chain
	for i in range(N):
		r1trial = r1 + h*np.random.uniform(0, 1) # Calculate a trial position 
		r2trial = r2 + h*np.random.uniform(0, 1)

		if metropolis(alpha, r1, r2, r1trial, r2trial) == True:
			r1 = r1trial # Keep the trial positions
			r2 = r2trial	
			# Accumulate the local energy, and local energy squared (for variance)
			E1 += localenergy(alpha, r1, r2)
			E2 += localenergy(alpha, r1, r2)**2
			steps += 1	

	# Divide by total number of steps to get energy
	tempenergy = E1/steps 
	tempvariance = tempenergy**2 - E2/steps

	return tempenergy, tempvariance

# Define a range of alpha to vary (Pade-Jastrow Variational Parameter)
alpharange = np.arange(0.2, 1.5, (1.5-0.2)/N)

# Initialize local energy and variance lists
energylist = []
variancelist = []

# Exact solution of the helium atom ground state from \cite{c4}
def exact(z):
	return z**2 - 4*z + 5/8*z

exactlist = []

# Calculate the ground state for varying alpha
for alpha in alpharange:
	energylist.append(main(alpha)[0])
	variancelist.append(main(alpha)[1])
	exactlist.append(exact(alpha))

if __name__ == "__main__": 
	plt.figure(1)
	plt.plot(alpharange, energylist, color = 'r', linestyle = '--', 
		label = 'VMC Simulation')
	plt.plot(alpharange, exactlist, label = 'Exact Solution')
	plt.xlabel(r'$\alpha$')
	plt.ylabel('Ground State Energy')
	plt.legend().draggable()

	plt.figure(2)
	plt.plot(alpharange, variancelist)
	plt.xlabel(r'$\alpha$')
	plt.ylabel('Variance of Energy')
	plt.show()