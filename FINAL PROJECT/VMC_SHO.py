#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Variational Monte Carlo Method for a 1D Simple
Harmonic Oscillator 
"""

__author__ = "Eric Yeung"

np.random.seed(1337) # Seed the random generator 

N = 1000 # Number of Monte-Carlo steps
h = 1    # Step-size

# Trial wave function
def psi(alpha, x):
	return np.exp(-2*alpha*x**2)

def localenergy(alpha, x):
	return alpha + (0.5 - 2*alpha**2)*x**2

# Metropolis Algorithm to keep position or accept a trial position 
def metropolis(alpha, x, xtrial):
	# Ratio between the probability functions |psi|^2
	omega = (psi(alpha, xtrial)/psi(alpha, x))**2 

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
	x = np.random.rand()

	# Initialize energy and energy squared
	E1 = 0.
	E2 = 0.

	# Keep track of steps
	steps = 0.

	# Begin the the loop for the Markov Chain
	for i in range(N):
		xtrial = x + h*np.random.uniform(0, 1) # Calculate a trial position 

		if metropolis(alpha, x, xtrial) == True:
			x = xtrial # Keep the trial position
			# Accumulate the local energy, and local energy squared (for variance)
			E1 += localenergy(alpha, x)
			E2 += localenergy(alpha, x)**2
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

# Calculate the ground state for varying alpha
for alpha in alpharange:
	energylist.append(main(alpha)[0])
	variancelist.append(main(alpha)[1])

if __name__ == "__main__": 
	plt.figure(1)
	plt.plot(alpharange, energylist, color = 'r', linestyle = '--', 
		label = 'VMC Simulation')
	plt.axhline(y = 0.5, label = "Exact solution")
	plt.xlabel(r'$\alpha$')
	plt.ylabel('Ground State Energy')
	plt.legend().draggable()

	plt.figure(2)
	plt.plot(alpharange, variancelist)
	plt.xlabel(r'$\alpha$')
	plt.ylabel('Variance of Energy')
	plt.show()
