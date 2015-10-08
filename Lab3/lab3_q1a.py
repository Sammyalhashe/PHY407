# PHY407, Fall 2015: Lab03, Q1a
# Author: Duong, Bang Chi

# Import modules
from gaussxw import gaussxw
from numpy import sin, cos, pi, sqrt, linspace, zeros
from matplotlib.pyplot import plot, title, xlabel, ylabel, grid, show

# Define constants
length_x = 70
x = linspace(-5,5,length_x) # in m
z = 3.0 # in m
wavelength = 1.0 # in m
N = 50 # number of sample points

# Define u
def u(x,z,wavelength):
	return x * sqrt(2.0/(wavelength*z))

# Define arrays for C(u) and S(u)
C_Array = zeros(length_x)
S_Array = zeros(length_x)

# Define C(u)
def C(u_func,x_i):
    # Define integrand of C(u)
    def f(t):
        return cos(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_i,z,wavelength)
	
    # Calculate the sample points and weights, then map them
    # to the required integration domain
    t, w = gaussxw(N)
    tp = 0.5*(t_final-t_initial)*t +0.5*(t_final+t_initial)
    wp = 0.5*(t_final-t_initial)*w

    # Do the integration
    integral_C = 0.0
    for k in range(N):
        integral_C += wp[k]*f(tp[k])
    return integral_C

# Define S(u)
def S(u_func,x_i):
    # Define integrand of S(u)
    def g(t):
        return sin(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_i,z,wavelength)

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    t, w = gaussxw(N)
    tp = 0.5*(t_final-t_initial)*t +0.5*(t_final+t_initial)
    wp = 0.5*(t_final-t_initial)*w
        
    # Do the integration
    integral_S = 0.0
    for k in range(N):
        integral_S += wp[k]*g(tp[k])
    return integral_S

# For each element in x, we have one C(u) and one S(u)
for i in range(length_x):
    C_Array[i] = C(u,x[i])
    S_Array[i] = S(u,x[i])

def I_by_I_zero(x):
	return (1.0/8.0) * ( (2*C_Array+1)**2 + (2*S_Array+1)**2 )

plot(x, I_by_I_zero(x))
title("I/I_0 vs x")
xlabel("x (m)")
ylabel("I/I_0")

show()
