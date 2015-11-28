# Use Hit/Miss Monte Carlo integration method to estimate integral of 1D
# function f(x) over specified range (a,b). f(x) must stay within (ymin,ymax)
# Oliver Watt-Meyer, November 16, 2015. Adapted from Newman solutions.
from math import sin,sqrt
from random import random

# define integrand
def f(x):
    return sin(1.0/(x*(2-x)))**2

N=10000 # number of points to use for Monte Carlo integration
count=0 # variable to keep track of points where y<f(x)

# domain of integration
a=0.0
b=2.0

# max and min values f(x) takes over range (a,b)
ymax=1.0
ymin=0.0

# 2D area
A=(b-a)*(ymax-ymin)

for i in range(N):
    # pick random point in plane (a,b) by (ymin,ymax)
    x=(b-a)*random()+a
    y=(ymax-ymin)*random()+ymin
    
    # check if point is below f(x) curve
    if y<f(x):
        count+=1

# compute integral and error using Equations 10.23 and 10.26
I=A*float(count)/N
sigma=sqrt(I*(A-I)/N)

print 'With Hit/Miss method, I = %.4f plus/minus %.4f'%(I,sigma)