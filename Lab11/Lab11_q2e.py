#!/usr/bin/env python
from __future__ import division
from Lab11_q2abc import *

"""
Ising Model, plot the up and down spins at 
T = 1, T = 2, T = 3, and an extreme case T = 300.
I don't have the visual module installed 
and plt.ion() is very buggy, so i'll just plot the final 
configuration of the dipoles
"""

__author__ = "Eric Yeung"

# Position of the spin up/down dipoles
xup, yup = np.where(spin == 1)      # Gives positions that yield +1 
xdown, ydown = np.where(spin == -1) # Gives positions that yield -1 

"""
plt.plot(xup, yup, linestyle = 'none', marker = r'$\uparrow$', 
	color = 'r', markersize = 15)
plt.plot(xdown, ydown, linestyle = 'none', marker = r'$\downarrow$', 
	color = 'b', markersize = 15)

Arrows, harder to see
"""

# From the question, I'll use squares. 
plt.plot(xup, yup, 'rs', markersize = 15, label = 'UP')
plt.plot(xdown, ydown, 'bs', markersize = 15, label = 'DOWN')

plt.title('Configuration of Dipoles in the Lattice at T = ' + str(temp))
plt.xlim(-1,dips) # Set the x resolution to include all 20 dipoles
plt.xticks([])
plt.ylim(-1,dips) # Set the y resolution to include all 20 dipoles 
plt.yticks([])

#plt.legend().draggable()
plt.show()