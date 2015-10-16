# PHY407, Fall 2015, Lab5 Q1b
# Author: DUONG, BANG CHI

#The scipy.io.wavfile allows you to read and write .wav files
from scipy.io.wavfile import read, write
from numpy import arange, linspace, empty, array
from numpy.fft import rfft, irfft
from matplotlib.pyplot import figure, subplot, plot, legend, xlabel, xlim, ylabel, show

#read the data into two stereo channels
#sample is the sampling rate, data is the data in each channel,
#  dimensions [2, nsamples]
sample, data = read('GraviteaTime.wav')
#sample is the sampling frequency, 44100 Hz
#separate into channels
channel_0 = data[:,0]
channel_1 = data[:,1]
N_Points = len(channel_0)
dt = 1/float(sample)
t0 = N_Points*dt
N1 = int(0.05/dt)
T = linspace(0.0,t0,N_Points)
T1 = linspace(0, 0.05, N1)

#... do work on the data...
fft_channel_0 = rfft(channel_0)
fft_channel_1 = rfft(channel_1)
freq = arange(N_Points/2+1)/t0
fc = 880.0

for i in range(len(fft_channel_0)):
    if freq[i]>fc:
        fft_channel_0[i] = 0
    else:
        fft_channel_0[i] = fft_channel_0[i]

for i in range(len(fft_channel_1)):
    if freq[i]>fc:
        fft_channel_1[i] = 0
    else:
        fft_channel_1[i] = fft_channel_1[i]
        
channel_0_out = irfft(fft_channel_0)
channel_1_out = irfft(fft_channel_1)

fft_channel_0_out = rfft(channel_0_out)
fft_channel_1_out = rfft(channel_1_out)

# Plot for input (whole time point)
figure(1)
subplot(211)
plot(T,array(channel_0),label='Original time series')
ylabel('Channel 0 data')
legend()

subplot(212)
plot(T,array(channel_1),label='Original time series')
legend()
xlabel('Time (seconds)')
ylabel('Channel 1 data')

# Plot for input vs filtered data (short time interval)
figure(2)
subplot(211)
plot(T1,array(channel_0[:N1]),label='Original time series')
plot(T1,array(channel_0_out[:N1]), 'r-', label='Filtered time series')
ylabel('Channel 0 data')
legend()

subplot(212)
plot(T1,array(channel_1[:N1]),label='Original time series')
plot(T1,array(channel_1_out[:N1]), 'r-', label='Filtered time series')
legend()
xlabel('Time (seconds)')
ylabel('Channel 1 data')

# Plot for Fourier coefficients
figure(3)
subplot(221)
plot(abs(fft_channel_0), label = 'Channel-0 Original Fourier coefficients')
legend()
xlim(0,8000)

subplot(222)
plot(abs(fft_channel_1), label = 'Channel-1 Original Fourier coefficients')
legend()
xlim(0,8000)

subplot(223)
plot(abs(fft_channel_0_out), label = 'Channel-0 Filtered Fourier coefficients')
legend()
xlim(0,8000)

subplot(224)
plot(abs(fft_channel_1_out), label = 'Channel-1 Filtered Fourier coefficients')
legend().draggable
xlim(0,8000)

#this creates an empty array data_out with the
#same shape as "data" (2 x N_Points) and the same
#type as "data" (int16)
data_out = empty(data.shape, dtype = data.dtype)
#fill data_out
data_out[:,0] = channel_0_out
data_out[:,1] = channel_1_out
write('GraviteaTime_lpf.wav', sample, data_out)
print "Finish filtering!"

show()
