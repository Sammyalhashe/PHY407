# PHY407, Fall 2015, Lab 10, Question 1
# Author: DUONG, BANG CHI

# Import modules
from __future__ import division, print_function
from pylab import *
from numpy.random import random
from visual import sphere

N = 500 # Number of points

# Generate random u and v
u = random(N)
v = random(N)

# Define theta and phi based on u and v
theta = arccos(1 - 2*u)
phi = 2*pi*v

# Define the coordinate of the points of the unit sphere
x = sin(theta)*cos(phi)
y = sin(theta)*sin(phi)
z = cos(theta)

position = zip(x,y,z)

# Run the visual module
for position_i in position:
	sphere(pos = position_i, radius=0.02)
