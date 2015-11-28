# PHY407, Fall 2015, Lab09, Question 2ab
# Author: DUONG, BANG CHI

# Import modules
from __future__ import division, print_function
from banded import banded
from pylab import *
from scipy.constants import hbar

# Define constants
h = 1e-18 # Size of time-step, in s
#hbar = 1.0546e-36 # Planck's constant, in J*s
L = 1e-8 # Length of the box, in m
M = 9.109e-31 # Mass of an electron, in kg
N = 1000 # Grid slices

a = L/N # Spacing of the spacial grid points


# Define components a1,a2 and b1,b2 of matrices A and B (both symmetric and tridiagonal) respectively
a1 = 1 + (h*1j*hbar / (2*M*a**2)) 
a2 = -h*1j*hbar / (4*M*a**2)
b1 = 1 - (h*1j*hbar / (2*M*a**2)) 
b2 =  h*1j*hbar / (4*M*a**2)

# Define a complex array for wavefunction psi 
psi = zeros(N+1, complex)

        # Define psi(time t=0, position x)
def psi0(x):
	x0 = L/2
	sigma = 1e-10 # in m 
	k = 5e10 # in m**-1
	return exp(-(x-x0)**2/2/sigma**2)*exp(1j*k*x)

# Define the position x
x = linspace(0, L, N+1)
psi[:] = psi0(x)
psi[[0,N]]=0

A = empty((3,N), complex)

A[0,:] = a2
A[1,:] = a1
A[2:,] = a2



#==============================================================================
# for i in range(100):
# 	v = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1])
# 	psi[1:N] = banded(A,v,1,1)
#==============================================================================

#plot(psi)

#==============================================================================
# Ap = zeros((N-1,N-1),complex)
# for i in range(N-2):
# 	Ap[i,i] = a2
# 	Ap[i+1,i] = a1 #Bottom
# 	Ap[i,i+1] = a1 #Right
# Ap[N-2,N-2] = a2
#==============================================================================

# Static plots at time t0
plot(real(psi), label = '$\Psi(x,t)$')
plot(abs(psi), label = '$|\Psi(x,t)|$')
plot(-abs(psi), label = '$-|\Psi(x,t)|$')
xlabel('$x$')
ylabel('$\Psi(x,t)$')
legend()
show()

# Animation
from visual import curve, rate

psi_c = curve()
psi_c.set_x(x-L/2)

#psi = banded(A,v,1,1)
while True:
	rate(30)
	psi_c.set_y(real(psi)*1e-9)
	psi_c.set_z(imag(psi)*1e-9)	
	for i in range(20):
		v = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1])
		psi[1:N] = banded(A,v,1,1)

