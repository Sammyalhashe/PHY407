# This script provides some plotting tips, using Question 2 from Lab 05. You 
# need to have SLP.txt, lon.txt and times.txt from Lab 05 in the same folder as
# the script.
# Oliver Watt-Meyer
# October 21, 2015

from numpy import loadtxt,arange,fft,copy
from matplotlib.cm import get_cmap
from pylab import figure,contourf,colorbar,xlabel,ylabel,title,savefig,show,\
     subplots_adjust,subplot,suptitle,tight_layout,yticks

# Load SLP and coordinates
SLP = loadtxt('SLP.txt')
Longitude = loadtxt('lon.txt')
Times = loadtxt('times.txt')

# Plot SLP
figure(1)
contourf(Longitude, Times, SLP)
xlabel('longitude(degrees)')
ylabel('days since Jan. 1 2015')
title('SLP anomaly (hPa)')
colorbar()

# Change figure size, yticks, colormap, contour intervals, add colorbar label and replot
my_cmap=get_cmap('RdBu_r')
my_contours=arange(-50,56,10)
time_ticks=[0,31,59,90]
time_labels=['Jan 1\n2015','Feb 1\n2015','Mar 1\n2015','Apr 1\n2015']

figure(2,figsize=(5,6))
contourf(Longitude, Times, SLP, my_contours, cmap=my_cmap)
cb=colorbar(ticks=my_contours)
cb.set_label('SLP anomaly (hPa)')
xlabel('Longitude (degrees)')
yticks(time_ticks,time_labels)
title('Sea level pressure anomaly at 50S')

# Compute rfft of SLP data over longitude axis to get m=3 and m=5 components
SLP_fft=fft.rfft(SLP)

for m in [3,5]:
    SLP_fft_m=copy(SLP_fft)
    SLP_fft_m[:,:m]=0.0
    SLP_fft_m[:,m+1:]=0.0
    
    SLP_m=fft.irfft(SLP_fft_m)
    
    if m==3:
        SLP_m3=copy(SLP_m)
    elif m==5:
        SLP_m5=copy(SLP_m)   
        
# define contour intervals for m=3 and m=5 plots
my_contours2=arange(-20,21,4)

# plot SLP data for all wavenumbers, m=3 and m=5
figure(3,figsize=(12,6))
suptitle('Sea level pressure anomaly at 50S (hPa)',fontsize=14)

subplot(131)
contourf(Longitude, Times, SLP, my_contours, cmap=my_cmap)
cb=colorbar(ticks=my_contours,use_gridspec=True)
xlabel('Longitude (degrees)')
yticks(time_ticks,time_labels)
title('All wavenumbers')

subplot(132)
contourf(Longitude, Times, SLP_m3, my_contours2, cmap=my_cmap)
cb=colorbar(ticks=my_contours2,use_gridspec=True)
xlabel('Longitude (degrees)')
yticks(time_ticks,time_labels)
title('m=3')

subplot(133)
contourf(Longitude, Times, SLP_m5, my_contours2, cmap=my_cmap)
cb=colorbar(ticks=my_contours2,use_gridspec=True)
xlabel('Longitude (degrees)')
yticks(time_ticks,time_labels)
title('m=5')

tight_layout() # note "use_gridspec=True" in colorbar() calls. Necessary for tight_layout to work with colorbars.
subplots_adjust(top=0.88) # adjust the top of the subplots to be lower so that suptitle fits nicely.

show()