# PHY407, Fall 2015, Lab6, Q1
# Author: DUONG, BANG CHI

# Import modules
from scipy.io.wavfile import read, write
from numpy import linspace, empty, pi
from matplotlib.pyplot import figure, subplot, plot, legend, xlabel, xlim, ylabel, show

#read the data into two stereo channels
#sample is the sampling rate, data is the data in each channel,
#  dimensions [2, nsamples]
sample, data = read('GraviteaTime.wav')
#sample is the sampling frequency, 44100 Hz
#separate into channels
channel_0_in = data[:,0]
channel_1_in = data[:,1]
data_length = len(channel_0_in)
dt = 1.0/float(sample)
t_final = data_length * dt
T = linspace(0, t_final, data_length)

N1 = int(0.05/dt)
T1 = linspace(0, 0.05, N1)

RC = (1.0/880.0) / 4.4 # in Hz

# Define the function f(V_out, V_in)
def f(channel_out_element, channel_in_element):
    return (channel_in_element - channel_out_element) / RC

channel_0_out = []
element_0_out = 0.0
for i in range(data_length):
    k1 = dt * f(element_0_out, channel_0_in[i])
    k2 = dt * f(element_0_out+0.5*k1, channel_0_in[i]+0.5*dt)
    k3 = dt * f(element_0_out+0.5*k2, channel_0_in[i]+0.5*dt)
    k4 = dt * f(element_0_out+k3, channel_0_in[i]+dt)
    element_0_out += (k1+2*k2+2*k3+k4)/6
    channel_0_out.append(element_0_out)

channel_1_out = []
element_1_out = 0.0
for i in range(data_length):
    k1 = dt * f(element_1_out, channel_1_in[i])
    k2 = dt * f(element_1_out+0.5*k1, channel_1_in[i]+0.5*dt)
    k3 = dt * f(element_1_out+0.5*k2, channel_1_in[i]+0.5*dt)
    k4 = dt * f(element_1_out+k3, channel_1_in[i]+dt)
    element_1_out += (k1+2*k2+2*k3+k4)/6
    channel_1_out.append(element_1_out)

# Plot for input (whole time point)
figure(1)
subplot(211)
plot(T, channel_0_in, label='Original time series')
plot(T, channel_0_out, label = 'Filtered time series') 
ylabel('Channel 0 data')
legend()

subplot(212)
plot(T, channel_1_in, label='Original time series')
plot(T, channel_1_out, label = 'Filtered time series') 
legend()
xlabel('Time (seconds)')
ylabel('Channel 1 data')

# Plot for input vs filtered data (short time interval)
figure(2)
subplot(211)
plot(T1, channel_0_in[:N1], label='Original time series')
plot(T1, channel_0_out[:N1], 'r-', label='Filtered time series')
ylabel('Channel 0 data')
legend()

subplot(212)
plot(T1, channel_1_in[:N1], label='Original time series')
plot(T1, channel_1_out[:N1], 'r-', label='Filtered time series')
legend()
xlabel('Time (seconds)')
ylabel('Channel 1 data')

	
# Create an empty array for data output	
data_out = empty(data.shape, dtype = data.dtype)
# Fill data_out
data_out[:,0] = channel_0_out
data_out[:,1] = channel_1_out
write('GraviteaTime_ODE.wav', sample, data_out)
print "Finished filtering, check the output file!"

show()
