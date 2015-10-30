# PHY407, Fall 2015: Lab03, Q1c
# Author: Duong, Bang Chi

# Import modules
from gaussxw import gaussxw
from numpy import sin, cos, pi, sqrt, linspace, zeros
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, grid, legend, contourf, show

# Define constants
length_x = 60
x = linspace(-5,5,length_x) # in m
z = linspace(1,5,length_x) # in m
wavelength = [1.0,2.0,4.0] # in m
N = 50 # number of sample points

# Define u
def u(x_value,z_value,wavelength_value):
	return x_value * sqrt(2.0/(wavelength_value*z_value))

# Define arrays for C(u) and S(u)
C_Array_1 = zeros((length_x,length_x))
S_Array_1 = zeros((length_x,length_x))

C_Array_2 = zeros((length_x,length_x))
S_Array_2 = zeros((length_x,length_x))

C_Array_4 = zeros((length_x,length_x))
S_Array_4 = zeros((length_x,length_x))

# Define C(u)
def C(u_func, x_value, z_value, wavelength_value, sample_points):
    # Define integrand of C(u)
    def f(t):
        return cos(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_value, z_value, wavelength_value)
	
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
def S(u_func, x_value, z_value, wavelength_value, sample_points):
    # Define integrand of S(u)
    def g(t):
        return sin(0.5 * pi * t**2)
    # Define boundaries of the integral, and number of sample points N
    t_initial = 0.0
    t_final = u_func(x_value, z_value, wavelength_value)

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
    for k in range(length_x):
        C_Array_1[i,k] = C(u, x[k], z[i], wavelength[0], N)
        S_Array_1[i,k] = S(u, x[k], z[i], wavelength[0], N)

        C_Array_2[i,k] = C(u, x[k], z[i], wavelength[1], N)
        S_Array_2[i,k] = S(u, x[k], z[i], wavelength[1], N)

        C_Array_4[i,k] = C(u, x[k], z[i], wavelength[2], N)
        S_Array_4[i,k] = S(u, x[k], z[i], wavelength[2], N)
		
# Define I/I_0 
def I_by_I_zero_1(x_value, z_value):
	return (1.0/8.0) * ( (2 * C_Array_1+1)**2 + (2 * S_Array_1+1)**2 )
	
def I_by_I_zero_2(x_value, z_value):
	return (1.0/8.0) * ( (2 * C_Array_2+1)**2 + (2 * S_Array_2+1)**2 ) 
	
def I_by_I_zero_4(x_value, z_value):
	return (1.0/8.0) * ( (2 * C_Array_4+1)**2 + (2 * S_Array_4+1)**2 ) 

# Plot the contours
figure(1)
contourf(z,x,I_by_I_zero_1(x, z), levels = linspace(0.0,1.4,15))
title("Contour of intensity when wavelength = 1.0 m")
xlabel("z (m)")
ylabel("x (m)")

figure(2)
contourf(z,x,I_by_I_zero_2(x, z), levels = linspace(0.0,1.4,15))
title("Contour of intensity when wavelength = 2.0 m")
xlabel("z (m)")
ylabel("x (m)")

figure(3)
contourf(z,x,I_by_I_zero_4(x, z), levels = linspace(0.0,1.4,15))
title("Contour of intensity when wavelength = 4.0 m")
xlabel("z (m)")
ylabel("x (m)")

show()
