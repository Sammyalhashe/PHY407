#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Ising Model: plots the total magnetization of the 
square lattice as a function of time for fixed 
temperatures. 
"""

__author__ = "Eric Yeung"

J = 1.0      # positive interaction constant
kb = 1.0     # Boltzmann constant
dips = 20.0  # Number of dipoles
N = int(1e6) # Number of monte carlo steps
Tc = 4*J/kb  # Curie Temperature (4 in our case) -> a lot of fluctuation and 0 mag
temp = Tc    # Temperature (1,2,3, 300) for next question

# Want an array to hold values of spin.
spin = np.random.random_integers(0, 1, size = (dips, dips)) # Initially all 0 or +1
spin[spin == 0] = -1 # Change the 0's to -1's as desired 

# Total energy of the entire system 
def totalenergy(spin):
	adjx = np.sum(spin[:-1]*spin[1:]) # Sum in one direction 
	# Take transpose of our spin array for the other direction
	adjy = np.sum(np.transpose(spin)[:-1]*np.transpose(spin)[1:]) 
	return -J*(adjx + adjy)

# Metropolis Algorithm to keep old energy or replace with new energy
def metropolis(E1, E2, temp):
	if E2 <= E1:
		return True # Keep new energy
	else:
		# Acceptance probability
		if np.random.rand() < np.exp(-(E2-E1)/(kb*temp)):
			return True # Keep new energy
		else:
			return False # Keep old energy

# Total magnetization of the system
def totalmagnet(spin):
	return np.sum(spin)

# Initialize energy and magnetization 
energy = [totalenergy(spin)]
magnet = [totalmagnet(spin)]

# Initialize tplot for plotting
tplot = [0]

# Start the loop for the Markov Chain
for i in range(N):
	# Random Coordinates of each dipole on the grid (0,0) - (20,20)
	x = np.random.randint(0, dips)
	y = np.random.randint(0, dips)
	
	Eold = totalenergy(spin)
	spin[x][y] *= -1 # Now flip a random value of spin, and calculate the new energy
	Enew = totalenergy(spin)
	print "Step| %s" % i
	if metropolis(Eold, Enew, temp) == True:
		energy.append(Enew)              # Accept new flipped spin
		magnet.append(totalmagnet(spin)) # Magnetization doesnt have a condition 
		tplot.append(i)                  # Keep track of time for each energy value
	else:
		spin[x][y] *= -1 # flip again for a new energy 

if __name__ == "__main__": 
	#plt.plot(tplot, energy)
	plt.plot(tplot, magnet)
	plt.xlabel('Time')
	plt.ylabel('Magnetization')
	plt.title('Total Magnetization of a square lattice using the Ising model')
	plt.show()
else:
	print "Imported stuff from previous question."