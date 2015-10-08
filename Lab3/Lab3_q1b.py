# PHY407, Fall 2015: Lab03, Q1a
# Author: Duong, Bang Chi

# Import modules
from gaussxw import gaussxw
from numpy import sin, cos, pi, sqrt, linspace, zeros
from matplotlib.pyplot import plot, title, xlabel, ylabel, grid, legend, show

# Define constants
length_x = 70
x = linspace(-5,5,length_x) # in m
z = 3.0 # in m
wavelength = 1 # in m
N_50 = 50 # number of sample points
N_100 = 100

# Define u
def u(x,z,wavelength):
	return x * sqrt(2.0/(wavelength*z))

# Define arrays for C(u) and S(u)
C_Array_N_50 = zeros(length_x)
S_Array_N_50 = zeros(length_x)
C_Array_2N_50 = zeros(length_x)
S_Array_2N_50 = zeros(length_x)

C_Array_N_100 = zeros(length_x)
S_Array_N_100 = zeros(length_x)
C_Array_2N_100 = zeros(length_x)
S_Array_2N_100 = zeros(length_x)

# Define C(u)
def C(u_func, x_i, sample_points):
    # Define integrand of C(u)
    def f(t):
        return cos(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_i,z,wavelength)
	
    # Calculate the sample points and weights, then map them
    # to the required integration domain
    t, w = gaussxw(sample_points)
    tp = 0.5*(t_final-t_initial)*t +0.5*(t_final+t_initial)
    wp = 0.5*(t_final-t_initial)*w

    # Do the integration
    integral_C = 0.0
    for k in range(sample_points):
        integral_C += wp[k]*f(tp[k])
    return integral_C

# Define S(u)
def S(u_func, x_i, sample_points):
    # Define integrand of S(u)
    def g(t):
        return sin(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_i,z,wavelength)

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    t, w = gaussxw(sample_points)
    tp = 0.5*(t_final-t_initial)*t +0.5*(t_final+t_initial)
    wp = 0.5*(t_final-t_initial)*w
        
    # Do the integration
    integral_S = 0.0
    for k in range(sample_points):
        integral_S += wp[k]*g(tp[k])
    return integral_S

# For each element in x, we have one C(u) and one S(u)
for i in range(length_x):
    C_Array_N_50[i] = C(u, x[i], N_50)
    S_Array_N_50[i] = S(u, x[i], N_50)
    C_Array_2N_50[i] = C(u, x[i], 2*N_50)
    S_Array_2N_50[i] = S(u, x[i], 2*N_50)

    C_Array_N_100[i] = C(u, x[i], N_100)
    S_Array_N_100[i] = S(u, x[i], N_100)
    C_Array_2N_100[i] = C(u, x[i], 2*N_100)
    S_Array_2N_100[i] = S(u, x[i], 2*N_100)

# Define I/I_0 for two cases
def I_by_I_zero_N_50(x_value):
	return (1.0/8.0) * ( (2 * C_Array_N_50+1)**2 + (2 * S_Array_N_50+1)**2 )

def I_by_I_zero_2N_50(x_value):
	return (1.0/8.0) * ( (2 * C_Array_2N_50+1)**2 + (2 * S_Array_2N_50+1)**2 )

def I_by_I_zero_N_100(x_value):
	return (1.0/8.0) * ( (2 * C_Array_N_100+1)**2 + (2 * S_Array_N_100+1)**2 )

def I_by_I_zero_2N_100(x_value):
	return (1.0/8.0) * ( (2 * C_Array_2N_100+1)**2 + (2 * S_Array_2N_100+1)**2 )

# Define error
error_N_50 = abs(I_by_I_zero_2N_50(x) - I_by_I_zero_N_50(x))

error_N_100 = abs(I_by_I_zero_2N_100(x) - I_by_I_zero_N_100(x))

plot(x, error_N_50, "-b", label = "N=50")
plot(x, error_N_100, "-r", label = "N=100")
legend(loc = "upper left")
title("error vs x")
xlabel("x (m)")
ylabel("error")

show()
