#!/usr/bin/env python
from __future__ import division
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

my_contours1=np.arange(-50,56,10)
my_contours2=np.arange(-20,21,4)

time_ticks=[0,31,59,90]
time_labels=['Jan 1\n2015','Feb 1\n2015','Mar 1\n2015','Apr 1\n2015']

if __name__ == "__main__":
    
    # Plots the first part (a)

    plt.figure(1)
    plt.contourf(Longitude, Times, SLP, my_contours1, cmap = 'RdBu_r') # The _r reverse the colors

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP anomaly (hPa)')
    plt.yticks(time_ticks,time_labels,rotation=45)
    plt.colorbar(ticks=my_contours1,use_gridspec=True)

    # Plots the second part (b)

    plt.figure(2)
    plt.contourf(Longitude, Times, extraction(3), my_contours2,  cmap = 'RdBu_r')

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP m = 3 (hPa)')
    plt.yticks(time_ticks,time_labels,rotation=45)
    plt.colorbar(ticks=my_contours2,use_gridspec=True)


    plt.figure(3)
    plt.contourf(Longitude, Times, extraction(5), my_contours2,  cmap = 'RdBu_r')

    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP m = 5 (hPa)')
    plt.yticks(time_ticks,time_labels,rotation=45)
    plt.colorbar(ticks=my_contours2,use_gridspec=True)


    plt.tight_layout()
    plt.show()

else:
    print "The plots are omitted."



