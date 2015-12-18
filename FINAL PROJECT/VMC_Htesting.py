#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Variational Monte Carlo Method for the Hydrogen Atom
"""

__author__ = "Eric Yeung"

np.random.seed(1337) # Seed the random generator 

N = 500 # Number of Monte Carlo steps
h = 1.00 # Step-size

# Trial wave function
def psi(alpha, r):
	return np.exp(-alpha*r)

# Local energy 
def localenergy(alpha, r):
	return -0.5*alpha**2 + alpha/r - 1/r

# Metropolis Algorithm to keep position or accept a trial position 
def metropolis(alpha, r, rtrial):
	# Ratio between the probability functions |psi|^2
	omega = (psi(alpha, rtrial)/psi(alpha, r))**2 

	if omega >= 1:
		return True # New position is accepted unconditionally
	else:
		if omega > np.random.rand(): # Acceptance probability
			return True # New position is accepted 
		else:
			return False # New position is not accepted, keep old 

# The main program
def main(alpha):
	# Radius, random number from 0 to 1
	r = np.random.rand()

	# Initialize energy and energy squared
	E1 = []
	E2 = []

	# Keep track of steps
	steps = 0.

	# Begin the the loop for the Markov Chain
	for i in range(N):
		rtrial = r + h*np.random.uniform(0, 1) # Calculate a trial position 

		if metropolis(alpha, r, rtrial) == True:
			r = rtrial # Keep the trial position
			# Accumulate the local energy, and local energy squared (for variance)
			E1.append(localenergy(alpha, r))
			E2.append(localenergy(alpha, r)**2)
			steps += 1

	# Divide by total number of steps to get energy
	tempenergy = np.mean(E1)
	tempvariance = tempenergy**2 - np.mean(E2)

	return tempenergy, tempvariance

# Define a range of alpha to vary (Pade-Jastrow Variational Parameter)
alpharange = np.arange(0.2, 1.0, (1.0-0.2)/N)

# Initialize local energy and variance lists
energylist = []
variancelist = []

# Calculate the ground state for varying alpha
for alpha in alpharange:
	energylist.append(main(alpha)[0])
	variancelist.append(main(alpha)[1])

if __name__ == "__main__": 
	plt.figure(1)
	plt.plot(alpharange, energylist, color = 'r', linestyle = '--', 
		label = 'VMC Simulation')
	plt.axhline(y = -0.5, label = "Exact solution")
	plt.ylabel('Ground State Energy')
	plt.legend().draggable()

	plt.figure(2)
	plt.plot(alpharange, variancelist)
	plt.xlabel('Alpha')
	plt.ylabel('Variance of Energy')
	plt.show()
