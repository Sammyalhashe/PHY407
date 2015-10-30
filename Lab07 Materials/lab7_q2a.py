# PHY407, Fall 2015, Lab7, Q2a
# Author: DUONG, BANG CHI

# Import modules
from numpy import linspace, zeros, sqrt, array
from matplotlib.pyplot import figure, plot, scatter, axhline, axvline, title, \
     xlabel, ylabel, grid, annotate, legend, show

# Define constants
G = 6.6738e-11 # Newton's G in m^3/(kg.s^2)
M = 1.9891e30 # mass of the Sun in kg
m = 5.9722e24 # mass of the Earth in kg

# Define an array for time, initial time, final time, time step
time_ini = 0.0 # in s
time_final = 3600*24*365*10 # 10 years in s
step = 3600 # 1 hour in s
time_length = int(((time_final - time_ini) / step) +1)

time = linspace(time_ini, time_final, time_length)


# Define arrays for position, velocity, and distance R
x = zeros(time_length); y = zeros(time_length);
vx = zeros(time_length); vy = zeros(time_length);

# Define initial conditions
x[0] = 1.4710e11; y[0] = 0;
vx[0] = 0; vy[0] = 3.0287e4;

# Verlet method
    # Half-step from intial conditions (outside the loop)
vx[0.5] = vx[0] + 0.5 * step * (-G*M*x[0]/sqrt(x[0]**2 + y[0]**2)**3)
vy[0.5] = vy[0] + 0.5 * step * (-G*M*y[0]/sqrt(x[0]**2 + y[0]**2)**3)
    # Loop 
for i in range(time_length - 1):

    x[i+1] = x[i] + step * vx[i+0.5]
    y[i+1] = y[i] + step * vy[i+0.5]

    kx = step * (-G*M*x[i+1]/sqrt(x[i+1]**2 + y[i+1]**2)**3)
    ky = step * (-G*M*y[i+1]/sqrt(x[i+1]**2 + y[i+1]**2)**3)
    
    vx[i+1] = vx[i+0.5] + 0.5 * kx
    vy[i+1] = vy[i+0.5] + 0.5 * ky
	
    vx[i+1.5] = vx[i+0.5] + kx
    vy[i+1.5] = vy[i+0.5] + ky

# Energy
    # KE in Joules
KE = 0.5 * m * (vx**2 + vy**2)
    # PE in Joules
PE = -G*M*m / sqrt(x**2 + y**2)
    # Total energy in Joules
Energy = KE + PE

# Plot
figure()
plot(x,y)
sun = array([0, 0]) # Place sun at origin
scatter(sun[0], sun[1], 400, c=["r"])
annotate('Sun', xy = (0,0))
scatter(x[0], y[0], 50, c=["y"])
annotate('Perihelion', xy = (x[0]+0.25e10,y[0]))
axhline(color = 'black')
axvline(color = 'black')
xlabel('x (m)')
ylabel('y (m)')
title('Orbit of the Earth')

figure()
plot(time,KE,label='Kinetic Energy')
plot(time,PE,label='Potential Energy')
plot(time,Energy,label='Total Energy')
legend().draggable()
xlabel('Time (s)')
ylabel('Energy (Joules)')
grid('on')

figure()
plot(time, Energy)
xlabel('Time (s)')
ylabel('Total Energy (Joules)')
grid('on')

show()
