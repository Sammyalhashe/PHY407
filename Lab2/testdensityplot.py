import numpy as np
import pylab as plt
from math import *

from scipy import special as sp

def I(r):
    if (r == 0):
        return 1/2   # Using the limit found in the textbook. 
    else:
        return (sp.jv(0,k*r)/(k*r))**2

wavelength = 500        # nano meters
k = 2*pi/wavelength

N = 1000     
stepsize = 0.95 # Increase resolution

Intensity = np.empty([N,N])

for i in range(N):
    y = stepsize*(i - N/2) # centre on image
    for j in range(N):
        x = stepsize*(j - N/2)
        r = np.sqrt((x)**2+(y)**2)
        Intensity[i,j] = I(r)

plt.imshow(Intensity, vmax= 0.08, cmap='hot')
plt.colorbar() 
plt.show()

