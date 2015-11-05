from numpy import empty,zeros,max, copy
from pylab import imshow,gray,show
from time import clock
# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
target = 1e-3   # Target accuracy (reduced from original problem)

jacobi = False
#using numpy slice to get vectorization
jacobi_numpy = True
# Create arrays to hold potential values
phi = zeros([M+1,M+1],float)
phi[0,:] = V
phiprime = empty([M+1,M+1],float)
startclock = clock()
# Main loop
delta = 1.0
while delta>target:
    if jacobi:
        # Calculate new values of the potential
        for i in range(M+1):
            for j in range(M+1):
                if i==0 or i==M or j==0 or j==M:
                    phiprime[i,j] = phi[i,j]
                else:
                    phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                                     + phi[i,j+1] + phi[i,j-1])/4
    if jacobi_numpy:
        #now do the same thing but replacing using numpy slice notation
        phiprime = copy(phi) #this copies the boundary conditions quickly
        #and this updates the interior points.
        phiprime[1:-1,1:-1]=(phiprime[2:,1:-1]+phiprime[:-2,1:-1]+phiprime[1:-1,2:]+phiprime[1:-1,:-2])/4

    
    # Calculate maximum difference from old values
    delta = max(abs(phi-phiprime))
    print 'delta: ', delta
            
    # Swap the two arrays around
    phi,phiprime = phiprime,phi    

endclock = clock()
time_clock = endclock - startclock
print 'time elapsed is %f seconds'%time_clock
# Make a plot
imshow(phi)
gray()
show()
