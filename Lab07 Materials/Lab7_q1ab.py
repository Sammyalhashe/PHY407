#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
Writes out the second order ODE in r 
as a pair of coupled First Order ODEs
"""

__author__ = "Eric Yeung"

from scipy.constants import epsilon_0, e, hbar, m_e

E_0 = 13.6 # Ionization energy of Hydrogen in eV
a_0 = 4*np.pi*epsilon_0*hbar**2/(m_e*e**2)
# Can decrease h and increase rmax to get more accurate soln
h = 0.002*a_0
rmax = 20*a_0 # Value at r -> infinity

def V(r):
	return -e**2/(4*np.pi*epsilon_0*r)

"""
I use a change in variable such that u = r*R, and the SODE becomes
u'' + (2*m/hbar**2*(E-V) - l*(l+1)/r**2)*u = 0
It's -now- trivial to reduce this to two first order ODEs
"""

def f(q, r, E): # General vector q
    u, v = q # Components of q
    du = v 
    dv = (l*(l+1)/r**2 - 2*m_e/hbar**2*(E-V(r)))*u 
    return np.array([du, dv], float)

def SchroSolve(E):
    """
    Initial Conditions at h (or very close to zero)
    Want psi(h) = 0, phi(h) = 1 or in our terms
    u(h) = 0, v(h) = 1
    """
    u = 0
    v = 1
    q = np.array([u, v], float)

    # Runge Kutta
    for r in np.arange(h, rmax, h):
        k1 = h*f(q, r, E)
        k2 = h*f(q + 0.5*k1, r + 0.5*h, E)
        k3 = h*f(q + 0.5*k2, r + 0.5*h, E)
        k4 = h*f(q + k3, r + h, E)
        q += (k1 + 2*k2 + 2*k3 + k4)/(6)

    # Recall the substitution R = u/r, where r = rmax is final step
    return q[0]/rmax

l = 1
n = 1
        
E_1 = -15*e/n**2
E_2 = -13*e/n**2
R2 = SchroSolve(E_1)

# Check with the analytical solution
theoryE1 = -13.6/1**2
theoryE2 = -13.6/2**2

accuracy = e/1000
while abs(E_1 - E_2) > accuracy:
    R1, R2 = R2, SchroSolve(E_2)
    E_1, E_2 = E_2, E_2 - R2*(E_2 - E_1)/(R2 - R1)

if __name__ == "__main__": 
        
    if n == 1:
        print "The ground state E_1 is %s eV for l = %s" % (E_2/e, l) # when n = 1
    
    elif n == 2:
        print "The first excited state E_2 is %s eV for l = %s" % (E_2/e, l) # when n = 2
    
    else:
        print "Invalid energy level, exiting."
        
else:
    print "Importing the energy levels and intial wavefunction equations"
