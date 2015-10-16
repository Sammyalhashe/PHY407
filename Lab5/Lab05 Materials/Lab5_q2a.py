#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

"""
Picks out the m =3, 5 fourier wavenumbers in SLP and plots them
"""

__author__ = "Eric Yeung"

SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

def extraction(m):

    finalSLP = np.empty(SLP.shape) 

    i = 0
    for i in range(len(Times)):
        temp = np.fft.rfft(SLP[i,:])
    
        for j in range(temp.shape[0]):
            if j != m:
                temp[j] = 0.           # Set each entry in each column(longitude) to 0 

        newSLP = np.fft.irfft(temp) 

        #Sets each row in finalSLP to the 1D array above 
    
        finalSLP[i,:] = newSLP 

    return finalSLP


if __name__ == "__main__":
    
    # Plots the first part (a)

    plt.figure(1)
    plt.contourf(Longitude, Times, SLP)

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP anomaly (hPa)')
    plt.colorbar()

    # Plots the second part (b)

    plt.figure(2)
    plt.contourf(Longitude, Times, extraction(3))

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP m = 3 (hPa)')
    plt.colorbar()


    plt.figure(3)
    plt.contourf(Longitude, Times, extraction(5))

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP m = 5 (hPa)')
    plt.colorbar()

    plt.tight_layout()
    plt.show()

else:
    print "The plots are omitted."



