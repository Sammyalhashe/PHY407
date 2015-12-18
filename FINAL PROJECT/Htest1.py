#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Variational Monte Carlo Method for the Hydrogen Atom
"""

__author__ = "Eric Yeung"

#np.random.seed(1337) # Seed the random generator 

N = 1000   # Number of Monte Carlo steps
h = 1       # Step-size

#alpha = 1. 



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
	e1 = 0.
	e2 = 0.

	# Keep track of steps
	steps = 0.

	# Begin the the loop for the Markov Chain
	for i in range(N):
		rtrial = r + h*np.random.uniform(0, 1) # Calculate a trial position 
		print rtrial
		
		if metropolis(alpha, r, rtrial) == True:
			r = rtrial # Keep the trial position
			# Accumulate the local energy, and local energy squared (for variance)
			e1 += localenergy(alpha, r)
			e2 += localenergy(alpha, r)**2
			steps += 1	

	tempenergy = e1/steps 
	tempvariance = tempenergy**2 - e2/steps

	return tempenergy, tempvariance

alpha = 0.2 # Pade-Jastrow Variational Parameter
print main(alpha)[0]

