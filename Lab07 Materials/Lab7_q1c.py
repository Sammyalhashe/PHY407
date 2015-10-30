#!/usr/bin/env python
from __future__ import division
from Lab7_q1ab import *

"""
Normalizes the eigenfunction R and plots it 
"""

__author__ = "Eric Yeung"

# Modified the program
def ModSolve(E):
    """
    Initial Conditions at h (or very close to zero)
    Want psi(h) = 0, phi(h) = 1 or in our terms
    u(h) = 0, v(h) = 1
    """
    u = 0
    v = 1
    q = np.array([u, v], float)
    Rwave = []

    # Runge Kutta
    for r in np.arange(h, rmax, h):
        k1 = h*f(q, r, E)
        k2 = h*f(q + 0.5*k1, r + 0.5*h, E)
        k3 = h*f(q + 0.5*k2, r + 0.5*h, E)
        k4 = h*f(q + k3, r + h, E)
        q += (k1 + 2*k2 + 2*k3 + k4)/(6)
        
        # Divide by 1e278 to get rid of some overflow errors and to get more points
        if np.isinf(q[0]/(rmax*1e278)) == True:
        	pass

        else:
        	Rwave.append(q[0]/(rmax*1e278))

    # Recall the substitution R = u/r, where r = rmax is final step
    return np.array([Rwave])

# l = 0, n =1 and l = 0, n = 2, and l = 1, n = 1
energys = np.array([-13.4999966366, -3.38780488905, -3.40127691697], float)

# Find our wavefunctions R(r) for each eigenenergy
R01 = ModSolve(energys[0])[0]
R02 = ModSolve(energys[1])[0]
R11 = ModSolve(energys[2])[0]

def simp(f):
	s = f[0] + f[-1] # First and final 

	for k in range(1, len(f), 2):
		s += 4*f[k+h]
	for k in range(2, len(f)-1, 2):
		s += 2*f[k*h]

	return s*h/3

# Normalization constant
unnorm01 = simp(abs(R01)**2); #norm01 = R01/unnorm01
unnorm02 = simp(abs(R02)**2); #norm02 = R02/unnorm02
unnorm11 = simp(abs(R11)**2); #norm11 = R11/unnorm11

# Apply normalization constant
norm01 = R01/np.sqrt(unnorm01)
norm02 = R02/np.sqrt(unnorm02)
norm11 = R11/np.sqrt(unnorm11)

rplot = np.arange(h, rmax, (rmax-h)/12)

if __name__ == "__main__": 

	plt.plot(rplot, norm01, label = 'Ground state l = 0')
	plt.plot(rplot, norm02, label = '1st Excited State l = 0')
	plt.plot(rplot, norm11, label = 'Ground state l = 1')
	#plt.xlim([0.8e-9, 1e-9])
	#plt.ylim([0, 10])
	
	plt.xlabel('r')
	plt.ylabel('R(r)')
	plt.title('Wavefunctions R(r) for different n, l')

	plt.legend().draggable()
	plt.show()

else:
	print "Importing the normalized wavefunctions R(r)"
