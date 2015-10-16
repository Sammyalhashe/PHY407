#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""
Calculates the energy eigenvalues of the quantum well
"""

__author__ = "Eric Yeung"

mmax = 100
nmax = 100
#a = 10 # electronvolts
a = 1.60218e-18 # joules 
L = 5e-10 # Meters
mass = 9.1094e-31 # kilograms
#hbar = 6.582119514e-16  # electronvolts * s
hbar = 1.054571800e-34  # joules * s

def Hcompute(m,n):

    if (m != n) and ((n%2 == 0 and m%2 == 0) or (n%2 != 0 and m%2 != 0)):
        Hvalue = 0
        #print "n and m either both odd or both even"

    elif (m != n) and ((n%2 == 0 and m%2 != 0) or (n%2 != 0 and m%2 == 0)):
        Hvalue = -1*(8*a*m*n)/(np.pi**2*(m**2 - n**2)**2)
        #print "n and m and odd/even in no order"

    elif (m == n):
        Hvalue = 0.5*a + (np.pi**2*hbar**2*m**2)/(2*mass*L**2)
        #print "n and m are equal"
    
    else:
        pass

    return Hvalue

def Hmatrix(m,n):
    Htemp = np.empty([mmax, nmax], float)
    
    for m in range(mmax):
        for n in range(nmax):
            Htemp[m,n] = Hcompute(m,n)
            
    return Htemp # This gives wrong ground state. 

def Hshift(m,n):
    Hshift = np.empty([mmax, nmax], float)
    
    for m in range(1, mmax + 1):
        for n in range(1, nmax + 1):
            Hshift[m-1,n-1] = Hcompute(m,n)
    
    return Hshift # symmetric and real eigenvalues = hermertian

if __name__ == '__main__':
	print "First %s Energy levels are:" % nmax 
	print np.linalg.eigvalsh(Hshift(0, 0))

	print "Ground state is given as %s, which is %s eV" % (np.linalg.eigvalsh(Hshift(0, 0))[0], 6.242e18*np.linalg.eigvalsh(Hshift(0, 0))[0])
	print "As expected. Done!"
    
else:
    print "_______Using the data in Lab4_q3c.py_______"
