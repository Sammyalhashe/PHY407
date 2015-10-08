# PHY407, Fall 2015: Lab03, Q1c
# Author: Duong, Bang Chi

# Import modules
from gaussxw import gaussxw
from numpy import sin, cos, pi, sqrt, linspace, zeros
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, colorbar, contourf, show

# Define constants
length_x = 100
x = linspace(-5,5,length_x) # in m
z = linspace(1,5,length_x) # in m
wavelength = [1.0,2.0,4.0] # in m
N = 100 # number of sample points

# Define arrays for C(u) and S(u)
C_Array = zeros((3,length_x,length_x))
S_Array = zeros((3,length_x,length_x))

# Define u
def u(x_value, z_value, wavelength_value):
    return x_value * sqrt(2.0/(wavelength_value * z_value))

# Define a general integrand
def Integrand(t, trig):
    return trig(0.5 * pi * t**2)

# Define a general trig function integral
def TrigIntegral(u_func, x_value, z_value, wavelength_value, sample_points, trig):
    # Define boundaries of the integral
    t_initial = 0.0
    t_final = u_func(x_value, z_value, wavelength_value)
	
    # Calculate the sample points and weights, then map them
    # to the required integration domain
    t, w = gaussxw(sample_points)
    wp = 0.5*(t_final-t_initial)*w
    tp = 0.5*(t_final-t_initial)*t +0.5*(t_final+t_initial)
    
    # Do the integration
    integral = 0.0
    for k in range(sample_points):
        integral += wp[k] * Integrand(tp[k], trig)
    return integral

# Define I/I_0 
def I_by_I_zero(x_value, z_value, wavelength_I):
	return (1.0/8.0) * ( (2.0 * C_Array[wavelength_I] + 1.0)**2.0 + (2.0 * S_Array[wavelength_I] + 1.0)**2.0 )

# Add values to the 3-d C and S arrays
for i in range(length_x):
    for k in range(length_x):
        for wavelengthI in range(3):
            C_Array[wavelengthI, i, k] = TrigIntegral(u, x[i], z[k], wavelength[wavelengthI], N, sin)
            S_Array[wavelengthI, i, k] = TrigIntegral(u, x[i], z[k], wavelength[wavelengthI], N, cos)
                    
# Plot the contours
for wavelengthI in range(3):
    figure(wavelengthI)
    contourf(z, x, I_by_I_zero(x, z, wavelengthI), levels = linspace(0.0,1.4,15))
    colorbar()
    title("Contour of intensity when wavelength = %s m" % (wavelength[wavelengthI]))
    xlabel("z (m)")
    ylabel("x (m)")
        
show()
