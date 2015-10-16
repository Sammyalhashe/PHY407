#various exercises with dft and pft
from numpy import zeros,loadtxt, arange, copy
from cmath import exp,pi
from pylab import plot,xlim,show, plot, title, subplot, figure, xlabel, savefig
from numpy.fft import rfft, irfft
from numpy import argsort
from time import time

#function to calculate the dft
def dft(y):
    N = len(y)
    c = zeros(N//2+1,complex)
    for k in range(N//2+1):
        for n in range(N):
            c[k] += y[n]*exp(-2j*pi*k*n/N)
    return c

#plot time series and dft
y = loadtxt("pitch.txt",float)
subplot(3,1,1)
plot(y)
title('pitch timeseries')
dft_time = time()
c = dft(y)
dft_time = time() - dft_time
subplot(3,1,2)
plot(abs(c))
title('amplitude of fourier coefficients')
xlim(0,500)

#------------------
#now do it again with FFT
fft_time = time()
c2 = rfft(y)
fft_time = time() - fft_time

#compare home made dft with fft performance
print 'dft time and fft time', dft_time, fft_time

#plot
subplot(3,1,3)
plot(abs(c2))
title('amplitude of fourier coefficients using fft')
xlim(0,500)

print 'maximum  |c2-c|: ', max(abs(c2-c))

#now do things with proper time dimensions and filter out desired frequencies
#sampling frequency for audio signal
f = 44100.0 #Hz
#related temporal sample
dt = 1/f #s
#length of vector
N = len(y)
#length of interval
T = N*dt
#convert to (angular) frequency
freq = arange(N/2+1)*2*pi/T
#dimensional time axis
t = arange(N)*dt
#sort on maximum frequency
MaxFreqs = argsort(abs(c2)**2) #get indexes of largest three frequencies
MaxFreqs = MaxFreqs[-1:-4:-1] #retain only top three
print 'top three frequencies and their amplitudes:'
print '%10.2f %10.2f %10.2f Hz'%(freq[MaxFreqs[0]]/(2*pi),freq[MaxFreqs[1]]/(2*pi),freq[MaxFreqs[2]]/(2*pi))
print '%10.2f %10.2f %10.2f'%(abs(c2[MaxFreqs[0]]),abs(c2[MaxFreqs[1]]), abs(c2[MaxFreqs[2]]))

#create a filtered array
c2_filt = copy(c2[:])
#zero out desired indices
c2_filt[MaxFreqs] = 0.0
#transform back to time domain
y_filt = irfft(c2_filt)
#now plot things dimensionally
figure(2,figsize=(6,10))
subplot(2,1,1)
plot(t,y, t, y_filt)
xlabel('t(s)')
title('pitch timeseries')
subplot(2,1,2)
plot(freq/(2*pi),abs(c2), freq/(2*pi), abs(c2_filt))
title('amplitude of fourier coefficients')
xlim((0,3000))
xlabel('f (Hz)')
show()
savefig('filtering_lab5.pdf')
figure(3)
#let's plot the cleaned up time series too - just for fun
plot(t,y-y_filt)
xlabel('t(s)')
title('pitch timeseries after removing lower amplitude signals')
