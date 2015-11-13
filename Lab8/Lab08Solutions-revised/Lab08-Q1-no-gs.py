#PHY407 2015, Lab 08, Q1
#From Newman Exercise 9.3
#import what's required from numpy
from numpy import zeros, copy, max, arange
from pylab import contour, show

#Set domain parameters
L = 0.1 #m
N = 100 
clo = N/5 #places the bottom of plate at 2 cm
chi = 4*N/5 #places the top of plate at 2 cm

#set up initial values and bioundary array
#set phi to zero everywhere
phi = zeros([N+1,N+1], float) #potential in volts
phi[clo:chi,clo] =  1.0 #V
phi[clo:chi,chi] = -1.0 #V

plates = zeros([N+1, N+1], int)
plates[clo:chi, clo] = 1#This creates a mask where the plates are located. These are taken as boundary
plates[clo:chi, chi] = 1 #conditions for the problems.

#solve for phi using Gauss-Seidel. Replace values on the fly.

delta = 1.0
target = 1e-5
iteration = 0
while delta>target:
    delta = 0.0
    oldphi = copy(phi)
    for i in range(1,N):
        for j in range(1,N):
            if plates[i,j] == 0:
                phi[i,j] = (oldphi[i+1,j]+oldphi[i-1,j]+oldphi[i, j+1]+oldphi[i,j-1])/4.0
    epsilon = abs(phi-oldphi).max()
    if epsilon>delta:
        delta = epsilon
        print 'for iteration %i delta is %f'%(iteration,delta)
        iteration+=1
contour(phi,colors='k',levels=arange(-1.0,1.1,0.1))
show()
                
