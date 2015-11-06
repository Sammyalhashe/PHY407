from __future__ import division, print_function

# PHY407, Fall 2015, Q3
__author__ = "DUONG, BANG CHI"

# Import modules
from numpy import zeros, sin, pi
from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, show, xlim


def T0(t):
    A = 10
    B = 12
    tau = 365 # days = 12 months
    return A + B*sin(2*pi*t/tau)

L = 20     # Depth in meters
D = 0.1   # Thermal diffusivity, m**2/days
N = 100       # Number of divisions in grid
a = L/N       # Grid spacing
h = 0.01     # Time-step
#epsilon = h/1000

# Define an array for temperature
T = zeros(N+1,float)

# Temperature everywhere equal to 10C, except at the surface and the deepest point
T[1:N]=10 

# Iteration
def iterate(T, t_min, t_max):
	# Main loop
	t = t_min
	c = h*D/a**2

	while t < t_max:
	
	    # Calculate the new values of T
		T[0] = T0(t)
		T[N] = 11
		T[1:N] = T[1:N] + c*(T[2:N+1]+T[0:N-1]-2*T[1:N])
	    
		t += h
	return T


T9 = iterate(T, 0, 365*9)

T9_i = T9
t_min = 365*9 # end of the 9th year

# Plot
for i in range(4):
	for t_max in [t_min + i*(365//4)]:
	#t_max = t_min + 365//4
		T9_i = iterate(T9_i, t_min, t_max)
		plot(T9_i, label = "the {}(st/nd/rd/th) 3-month intervals of the 10th year".format(i+1))
		t_min = t_max

legend(prop={'size':11}).draggable()
xlabel("Depth (m)")
ylabel("Temperature (degree C)")
xlim(0,20)
title("Thermal Diffusion in Earth's Crust")

show()
