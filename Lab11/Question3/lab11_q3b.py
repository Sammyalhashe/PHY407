# PHY407, Fall 2015, Lab 11, Question 3
# Author: DUONG, BANG CHI

# Import modules
from __future__ import division, print_function
from pylab import *
from numpy.random import random,randint,seed

# Define constants
Tmax = 10
Tmin = 1e-10
tau = 1e4
sigma = 1

# Define initial points
x0 = 2
y0 = 2

# Define lists for x,y, and t
t_list = []
x_list = []
y_list = []

n = randint(1,1e6)
seed(665600); print('Seed n = {0}'.format(n))
# Function to generate two Gaussian random numbers
def gaussian():
    r = sqrt(-2 * sigma**2 * log(1-random()))
    theta = 2*pi*random()
    x = r*cos(theta)
    y = r*sin(theta)
    return x,y

# Define functions g(x,y), f(x,y)
def g(x,y):
    if 0 < x < 50 and -20 < y < 20:
        return cos(x) + cos(sqrt(2)*x + cos(sqrt(3)*x) + (y-1)**2)
    else:
        return 1e10
	
def swap_function(f):	
	
    def h(x,y):
        return g(x,y)
	
    return h
	
'''Turn the @swap_function 'on' to get g(x,y), and 'off' to get f(x,y)'''
@swap_function
def f(x,y):
    return x**2 - cos(4*pi*x) + (y-1)**2
'''Turn the @swap_function 'on' to get g(x,y), and 'off' to get f(x,y)'''

function = f(x0,y0)
t = 0
T = Tmax
x, y = x0, y0

while T > Tmin:
	
    t+=1
    T = Tmax*exp(-t/tau)
	
    old_x, old_y = x, y
    old_function = function
    
    delta_x, delta_y = gaussian()
    x += delta_x
    y += delta_y
    function = f(x,y)

    delta_function = function - old_function
    
    if random() > exp(-delta_function/T):
        x, y = old_x, old_y
        function = old_function

    x_list.append(x)
    y_list.append(y)
    t_list.append(t)

print('(x,y) = ({},{}) with f(x,y) = {}'.format(x,y,function))

# Plots
figure()
scatter(t_list, x_list)
xlabel('time point')
ylabel('x coordinate')

figure()
scatter(t_list, y_list)
xlabel('time point')
ylabel('y coordinate')

show()
