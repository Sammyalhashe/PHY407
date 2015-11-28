from __future__ import division
from random import random,randrange
import numpy as np
import matplotlib.pyplot as plt

"""
Adapted from mcsim.py and the code given in the lab.
"""

__author__ = "Eric Yeung"

N = 1000
steps = int(2.5e5)

def testing(T):
	# Create a 2D array to store the quantum numbers
	n = np.ones([N,3],int)

	# Main loop
	E = 3*N*np.pi*np.pi/2
	for k in range(steps):

	    # Choose the particle and the move
	    i = randrange(N)
	    j = randrange(3)
	    if random()<0.5:
	        dn = 1
	        dE = (2*n[i,j]+1)*np.pi*np.pi/2
	    else:
	        dn = -1
	        dE = (-2*n[i,j]+1)*np.pi*np.pi/2

	    # Decide whether to accept the move
	    if n[i,j]>1 or dn==1:
	        if random() < np.exp(-dE/T):
	            n[i,j] += dn
	            E += dE

	#This calculates the energy of each particle,
	# neglecting constant factors
	energy_n = n[:,0]**2+n[:,1]**2+n[:,2]**2
	#This calculates the frequency distribution and creates a plot
	plt.figure(1)
	plt.clf()
	hist_output = plt.hist(energy_n,50)
	#This is the frequency distribution
	energy_frequency = hist_output[0]
	#This is what the x-axis of the plot should look like
	# if we plot the energy distribution as a function of n
	energy_vals = 0.5*(hist_output[1][:-1]+hist_output[1][1:])
	n_vals = np.sqrt(energy_vals)
	#Create the desired plot
	plt.figure(2)
	plt.clf()
	plt.bar(n_vals,energy_frequency, width = 0.1)
	plt.close(1); plt.close(2)

	navg = np.sum(hist_output[0]*n_vals)/np.sum(hist_output[0])

	return E, navg

#temprange = [10.0, 40.0, 100.0, 400.0, 1200.0, 1600.0]
temprange = np.arange(10.0, 1600.0, 50)

energys = []
averages = []

for T in temprange:
	energys.append(testing(T)[0])
	averages.append(testing(T)[1])
	print T, testing(T)[0], testing(T)[1] 

# Estimate the heat capacity 
Cv = (energys[-1] - energys[0])/(temprange[-1] - temprange[0])
print "Over this range, heat capacity is %.2f" % Cv

graph, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(temprange, energys, label = 'Stuff1')
ax1.set_xlabel('T')
ax1.set_ylabel('E(T)')

ax2.plot(temprange, averages, label = 'Stuff2')
ax2.set_xlabel('T')
ax2.set_ylabel(r'$\overline{n}$(T)')

plt.tight_layout()
plt.show()
