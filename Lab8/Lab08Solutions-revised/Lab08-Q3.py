#PHY407 2015, Lab 08, Q3
#From Newman Exercise 9.4
#import what's required from numpy
from numpy import empty, linspace, pi, sin
from pylab import plot, show, xlabel, ylabel

#Define constants and numerical parameters.
L = 20.0 #length of system in m
N = 100 #number of grid cells
a = L/N #Point separation
h = 0.01 #time spacing in days
D = 0.1 # Thermal diffusivity of crust in m^2/day

#constants for sinusoidal variation
A = 10.0
B = 12.0
tau = 365
Tfixed = 11.0 #Temperature at bottom

tmax = 10.01*365
steps = int(tmax/h)

#Just assume three month intervals are 0.25 years
s1 = int(9.25*365/h)
s2 = int(9.50*365/h)
s3 = int(9.75*365/h)
s4 = int(10.00*365/h)

#initialize T
T = empty(N+1, float)
T[0:N] = A
T[N] = Tfixed

#main loop
x = linspace(0.0, L, N+1)
for k in range(steps):
    t = k*h
    T[0] = A + B*sin(2*pi*t/tau)
    T[1:N] += h*D*(T[0:N-1]+T[2:N+1]-2*T[1:N])/(a**2)
    if k==s1 or k ==s2 or k==s3 or k==s4:
        plot(x,T)

xlabel("Depth(meters)")
ylabel("Temperature (Celsius)")
show()
