#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time

"""
Simulates linear water waves at various different times
using the spectral method 
"""

__author__ = "Eric Yeung"

L = 400.
D = 50. 
bdelta = 2.
A = 1.
g = 9.8 

h = 0.1
dx = 1.0
dz = 1.0

xpoints = np.arange(0, L, dx)
zpoints = np.arange(-D, 0, dz)

# Start timing 
startTime1 = time.time()

# Initialize eta
eta = np.zeros(len(xpoints), float)

# Initialize phi
phi = np.zeros([len(zpoints), len(xpoints)], float)

# IC for eta
for i in range(len(xpoints)):
    eta[i] = -A*np.exp(-(xpoints[i]-L/2)**2/bdelta**2)

# Create arrays to store our solutions in time
etavalues = np.zeros([len(xpoints), L], float)
etavalues [0,:] = eta

phivalues = np.zeros([len(xpoints), L], float)
phivalues [0,:] = -g*eta # append initial conditions at z = 0

nsteps = 200
t = 0
i = 1
# Time Evolution
while i < nsteps:
    """
    Don't care about k = 0 solution, because
    it does not evolve in time. Start k at 1.
    """
    karray = np.arange(1, len(xpoints)/2+2)*2*np.pi/L 
    omegak = np.sqrt(g*karray)

    etahat = np.fft.rfft(eta)
    etahat = etahat*karray*np.cos(omegak*t)
    etaprime = np.fft.irfft(etahat)
    etavalues[i,:] = etaprime

    phihat = np.fft.rfft(-g*eta)
    phihat = phihat*karray*np.sin(omegak*t)*np.exp(karray*0)/omegak # Iuno about this
    phiprime = np.fft.irfft(phihat)
    phivalues[i,:] = phiprime

    # Find the indices for different times
    #print t, i
    t += 2*h
    i += 1

etaplot0 = etavalues[0,:]
etaplot2 = etavalues[10,:]
etaplot10 = etavalues[50,:]
etaplot40 = etavalues[199,:]

phiplot0 = phivalues[0,:]
phiplot2 = phivalues[10,:]
phiplot10 = phivalues[50,:]
phiplot40 = phivalues[199,:]

if __name__ == "__main__":
    figureeta, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

    ax1.plot(xpoints, etaplot0)
    ax2.plot(xpoints, etaplot2)
    ax3.plot(xpoints, etaplot10)
    ax4.plot(xpoints, etaplot40)
    
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlabel('x')
        ax.set_ylabel('$\eta$(x,t)')

    ax1.set_title('Time = 0s')
    ax2.set_title('Time = 2s')
    ax3.set_title('Time = 10s')
    ax4.set_title('Time = 40s')

    plt.tight_layout()

    figurephi, ((ax5, ax6), (ax7, ax8)) = plt.subplots(2,2)

    ax5.plot(xpoints, phiplot0, color = 'r')
    ax6.plot(xpoints, phiplot2, color = 'r')
    ax7.plot(xpoints, phiplot10, color = 'r')
    ax8.plot(xpoints, phiplot40, color = 'r')
    
    """
    Contour Plot?

    ax5.contourf(xpoints, zpoints, phiplot0)
    ax6.contourf(xpoints, zpoints, phiplot2)
    ax7.contourf(xpoints, zpoints, phiplot10)
    ax8.contourf(xpoints, zpoints, phiplot40)
    """

    for ax in (ax5, ax6, ax7, ax8):
        ax.set_xlabel('x')
        ax.set_ylabel('$\phi$(x, z, t)')    

    ax5.set_title('Time = 0s')
    ax6.set_title('Time = 2s')
    ax7.set_title('Time = 10s')
    ax8.set_title('Time = 40s')
    
    time1 = time.time() - startTime1 # Stop timing second operation
    print "This simulation (spectral method) took %s seconds!" % time1
    
    plt.tight_layout()
    plt.show()
