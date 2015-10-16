#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from Lab4_q3c import *

"""
Plots the wavefunction and normalizes it. 
"""

__author__ = "Eric Yeung"

# Find our wavefunction from the eigenvectors for a given energylevel

E_n, psi_n = np.linalg.eigh(Hshift(0, 0))

def psi(energylevel, x):
    wavefunc = 0.
    for n in range(nmax):
        wavefunc += psi_n[energylevel,n]*np.sin(n*np.pi*x/L)
    return wavefunc

def Modulus_Squared(energylevel, x):
	return np.ma.conjugate(psi(energylevel, x))*psi(energylevel, x)

# Use simpson's rule for simplicity

x0 = 0.0
xf = L
h = (xf - x0)/nmax
energylevel = 1

simp = Modulus_Squared(energylevel, x0) + Modulus_Squared(energylevel, xf)

for k in range(1,nmax,2):
    simp += 4*Modulus_Squared(energylevel, x0+k*h)

for k in range(2,nmax-1,2):
    simp += 2*Modulus_Squared(energylevel, x0+k*h)
 
Unnormalized = simp*h/3; print "Unnormalized integral is %s" % Unnormalized

NormalizationConstant = np.sqrt(Unnormalized)
print "Normalization constant is %s"  % NormalizationConstant

def NormModSquared(energylevel, x):
	return Modulus_Squared(energylevel, x)/NormalizationConstant

# Plot the probability density of the ground state

xplot = np.linspace(0, L, len(psi_n))
plt.plot(xplot,NormModSquared(1, xplot), color = 'deepskyblue', label = 'Ground State')
plt.plot(xplot,NormModSquared(2, xplot), color = 'crimson', label = '1st Excited State')
plt.plot(xplot,NormModSquared(3, xplot), color = 'darkgreen', label = '2nd Excited State')

plt.xlabel('x')
plt.ylabel('$|\psi(x)|^2$')
plt.title('Probability Density $|\psi(x)|^2$')
plt.xlim(0,5e-10)


plt.legend().draggable()
plt.show()
